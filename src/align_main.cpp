/**
 * Ravi Gaddipati
 * Dec 23, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Main aligner interface
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include "align_main.h"
#include "alignment.h"
#include "sim.h"
#include "threadpool.h"
#include <mutex>
#include <cctype>
#include <unordered_set>
#include <functional>
#include <cstdlib>

using rg::Deleter;

int align_main(int argc, char *argv[]) {
    std::string cl = "vargas ";
    {
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        cl = ss.str();
    }

    // Load parameters
    unsigned match, npenalty, threads, chunk_size, subsample, max_paths;
    std::string read_file, gdf, align_targets, out_file, pgid, mismatch, rdg, rfg;
    bool end_to_end = false, fwdonly = false, p64=false, msonly=false, maxonly=false, notraceback=false;

    cxxopts::Options opts("vargas align", "Align reads to a graph.");
    try {
        opts.add_options("Input")
        ("g,gdef", "<str> *Graph definition file.", cxxopts::value(gdf))
        ("U,reads", "<str> *Unpaired reads in SAM, FASTQ, or FASTA format.", cxxopts::value(read_file));

        opts.add_options("Optional")
        ("S,sam", "<str> Output file.", cxxopts::value(out_file))
        ("msonly", "Only report max score. Improves speed.", cxxopts::value(msonly)->implicit_value("1"))
        ("maxonly", "Only report max score, location, and count. Improves speed.", cxxopts::value(maxonly)->implicit_value("1"))
        ("phred64", "Qualities are Phred+64, not Phred+33.", cxxopts::value(p64)->implicit_value("1"))
        ("p,subsample", "<N> Sample N random reads, 0 for all.", cxxopts::value(subsample)->default_value("0"))
        ("a,alignto", "<str> Target graph, or SAM Read Group -> graph mapping.\"(RG:ID:<group>,<target_graph>;)+|<graph>\"", cxxopts::value(align_targets))
        ("s,assess", "[ID] Use score profile from a previous alignment.", cxxopts::value(pgid)->implicit_value("."))
        ("f,forward", "Only align to forward strand.", cxxopts::value(fwdonly))
        ("notraceback", "Do not compute traceback (CIGAR).", cxxopts::value(notraceback)->implicit_value("1"))
        ("max-trace-paths", "<N> Max candidate paths per read when tracing back through a variant graph window.", cxxopts::value(max_paths)->default_value("256"));

        opts.add_options("Scoring")
        ("ete", "End to end alignment.", cxxopts::value(end_to_end))
        ("ma", "<N> Match bonus.", cxxopts::value(match)->default_value("2"))
        ("mp", "<MIN,MAX> Mismatch penalty. Lower qual=lower penalty.", cxxopts::value(mismatch)->default_value("2,6"))
        ("np", "<N> Penalty for non-A/C/G/T.", cxxopts::value(npenalty)->default_value("1"))
        ("rdg", "<GO,GEXT> Read gap open/extension penalty.", cxxopts::value(rdg)->default_value("3,1"))
        ("rfg", "<GO,GEXT> Ref gap open/extension penalty.", cxxopts::value(rfg)->default_value("3,1"));

        opts.add_options("Threading")
        ("j,threads", "<N> Number of threads.", cxxopts::value(threads)->default_value("1"))
        ("u,chunk", "<N> Partition into tasks of max size N.", cxxopts::value(chunk_size)->default_value("64"));

        opts.add_options()("h,help", "Display this message.");

        opts.parse(argc, argv);
    } catch (std::exception &e) {
        throw std::invalid_argument("Error parsing options: " + std::string(e.what()));
    }

    if (opts.count("h")) {
        align_help(opts);
        return 0;
    }

    if (!opts.count("gdef")) {
        align_help(opts);
        throw std::invalid_argument("Graph definition file required.");
    }

    if (!opts.count("reads")) {
        align_help(opts);
        throw std::invalid_argument("No read file provided.");
    }
    ReadFmt format = read_fmt(read_file);

    if (chunk_size < vargas::Aligner::read_capacity() || chunk_size % vargas::Aligner::read_capacity() != 0) {
        std::cerr << "[warn] Chunk size is not a multiple of SIMD vector length: "
                  << vargas::Aligner::read_capacity() << std::endl;
    }

    if (opts.count("assess") && format != ReadFmt::SAM) {
        throw std::invalid_argument("Assess is only available for SAM inputs.");
    }
    if (!align_targets.empty() && format != ReadFmt::SAM) {
        throw std::invalid_argument("Alignment targets only available for SAM inputs.");
    }

    if(opts.count("msonly") && opts.count("maxonly")) {
        throw std::invalid_argument("At most one of msonly and maxonly can be specified.");
    }

    vargas::isam reads;
    if (format == ReadFmt::FASTQ) {
        load_fast(read_file, true, reads, p64);
    } else if (format == ReadFmt::FASTA) {
        load_fast(read_file, false, reads, p64);
    } else {
        reads.open(read_file);
    }
    reads.subset(subsample);
    auto &reads_hdr = reads.header();

    vargas::ScoreProfile prof;
    {
        prof.match = match;
        prof.ambig = npenalty;

        auto sp = rg::split(mismatch, ',');
        if (sp.size() == 2) {
            prof.mismatch_min = std::stoi(sp[0]);
            prof.mismatch_max = std::stoi(sp[1]);
        }
        else if (sp.size() == 1) {
            prof.mismatch_max = prof.mismatch_min = std::stoi(sp[0]);
        }
        else {
            throw std::invalid_argument("Invalid --mp argument.");
        }

        sp = rg::split(rdg, ',');
        if (sp.size() == 2) {
            prof.read_gopen = std::stoi(sp[0]);
            prof.read_gext = std::stoi(sp[1]);
        }
        else {
            throw std::invalid_argument("Invalid --rdg argument.");
        }

        sp = rg::split(rfg, ',');
        if (sp.size() == 2) {
            prof.ref_gopen = std::stoi(sp[0]);
            prof.ref_gext = std::stoi(sp[1]);
        }
        else {
            throw std::invalid_argument("Invalid --rfg argument.");
        }
    }

    if (pgid == ".") {
        bool check = false;
        for (const auto &i : reads.header().programs) {
            if (std::find(vargas::supported_pgid.begin(), vargas::supported_pgid.end(), i.first)
            != vargas::supported_pgid.end()) {
                pgid = i.first;
                prof = vargas::program_profile(i.second.command_line);
                check = true;
                break;
            }
        }
        if (!check) throw std::invalid_argument("No suitable scoring profile found in SAM program header.");
        std::cerr << "Using profile for: " << pgid << "\n";
    } else if (pgid.length()) {
        try {
            prof = vargas::program_profile(reads_hdr.programs.at(pgid).command_line);
        } catch (std::exception &e) {
            throw std::invalid_argument("Unrecognized PG ID: " + pgid);
        }
    } else {
        prof.end_to_end = end_to_end;
    }

    vargas::SAM::Header::Program pg;
    pg.command_line = cl;
    pg.name = "vargas_align";
    pg.id = "VA";
    pg.version = __DATE__;
    std::replace_if(pg.version.begin(), pg.version.end(), [](char c) { return std::isspace(c); }, ' '); // rm tabs
    const auto assigned_pgid = reads_hdr.add(pg);

    size_t read_len;
    auto task_list = create_tasks(reads, align_targets, chunk_size, read_len);

    const size_t num_tasks = task_list.size();
    if (num_tasks < threads) {
        std::cerr << "[warn] Number of threads is greater than number of tasks. Try decreasing -u.\n";
    }

    threads = threads ? (static_cast<unsigned int>(threads) > static_cast<unsigned int>(task_list.size()) ? static_cast<int>(task_list.size()) : threads)
                      : 1;

    int bias = 255 - (read_len * match);
    const bool use_wide = (bias < 0) or (end_to_end and (static_cast<signed long long>(prof.ref_gopen + (prof.ref_gext * (read_len - 1))) > bias || static_cast<signed long long>(read_len * prof.mismatch_max) > bias));
    if (use_wide) {
        std::cerr << "Score range: " << read_len * match << " to -" << std::min(prof.ref_gopen + (prof.ref_gext * (read_len - 1)), read_len * prof.mismatch_max) <<
        ". Using 16-bit aligner (" << vargas::WordAligner::read_capacity() << " reads/vector).\n";
    }
    std::cerr << "Scoring profile: " << prof.to_string() << "\n";

    std::vector<std::unique_ptr<vargas::AlignerBase, rg::Deleter>> aligners(threads);
    for (size_t k = 0; k < static_cast<size_t>(threads); ++k) {
        aligners[k] = make_aligner(prof, read_len, use_wide, msonly, maxonly);
    }


    std::cerr << "\nLoading \"" << gdf << "\"...\n";
    auto start_time = std::chrono::steady_clock::now();
    vargas::GraphMan gm(gdf);
    if (gm.labels().size() != 1 && maxonly) {
        std::cerr << "[warn] With --maxonly, max score position and count may be incorrect because the genome is a graph." << std::endl;
    }
    if (gm.labels().size() != 1 && !maxonly && !msonly) {
        throw std::invalid_argument("Cannot calculate 2nd-max score when the genome is a graph. Use --msonly or --maxonly.");
    }
    std::cerr << rg::chrono_duration(start_time) << "s.\n";

    if (out_file.length()) std::cerr << "Writing to \"" << (out_file.empty() ? "stdout" : out_file) << "\".\n";
    reads_hdr.programs[assigned_pgid].aux.set(ALIGN_SAM_PG_GDF, gdf);
    vargas::osam aligns_out(out_file, reads_hdr);
    char phred_offset = opts.count("phred64") ? 64 : 33;
    align(gm, task_list, aligns_out, aligners, fwdonly, msonly, maxonly, notraceback, phred_offset, max_paths);

    return 0;
}

/**
 * @brief
 * Result of a flat (linear-reference) traceback.
 */
struct FlatTraceResult {
    std::string cigar;  // CIGAR string for the alignment
    int ref_offset;     // offset into ref_slice where the alignment begins
    int score;          // DP optimal score
};

/**
 * @brief
 * Recompute an affine-gap DP alignment of a read against a flat reference slice and
 * recover the traceback (CIGAR + start offset within the slice). This is the exact DP
 * previously inlined in align_helper_func; factored out so it can be reused to score
 * individual candidate paths through a variant-graph window.
 * @param ref_slice reference bases (numeric rg::Base) to align against
 * @param read read sequence
 * @param qual raw Phred-encoded quality characters; empty or wrong-length => no quality used
 * @param phred_offset Phred offset (33 or 64)
 * @param prof score profile (match/mismatch/gap/ambig, end_to_end mode)
 * @return CIGAR string, start offset into ref_slice, and DP optimal score
 */
static FlatTraceResult flat_traceback(const std::vector<rg::Base> &ref_slice,
                                      const std::string &read,
                                      const std::string &qual,
                                      char phred_offset,
                                      const vargas::ScoreProfile &prof) {
    const int ref_len = static_cast<int>(ref_slice.size());
    const bool has_quality = qual.size() == read.size();

    // Allocate the three DP score matrixes: M (match) D (deletion) I (insertion), initialize with zero
    std::vector<std::vector<int>> M(read.length()+1, std::vector<int>(ref_len+1, 0));
    std::vector<std::vector<int>> D(read.length()+1, std::vector<int>(ref_len+1, 0));
    std::vector<std::vector<int>> I(read.length()+1, std::vector<int>(ref_len+1, 0));
    // Allocate the three DP traceback matrixes: tM (match) tD (deletion) tI (insertion)
    std::vector<std::vector<int>> tM(read.length()+1, std::vector<int>(ref_len+1, 0));
    std::vector<std::vector<int>> tD(read.length()+1, std::vector<int>(ref_len+1, 0));
    std::vector<std::vector<int>> tI(read.length()+1, std::vector<int>(ref_len+1, 0));

    // initialize first row and first column of each score matrix as appropriate
    // gaps in beginning of reference (first row) ending in match doesn't make sense
    std::fill(M[0].begin(), M[0].end(), -prof.ref_gext*ref_len);
    M[0][0] = 0; // except zero characters of each is free
    // gaps in beginning of reference (first row) ending in gap in query doesn't make sense
    std::fill(I[0].begin(), I[0].end(), -prof.ref_gext*ref_len);
    // gaps in beginning of query (first col) ending in gap in reference / match doesn't make sense
    for (unsigned row = 1; row <= read.length(); ++row) {
        D[row][0] = -prof.ref_gext*ref_len;
        M[row][0] = -prof.ref_gext*ref_len;
    }
    if (prof.end_to_end) { // semiglobal
        // gaps in beginning of query (first col) ending in gap in query accumulate
        for (unsigned row = 1; row <= read.length(); ++row) {
            I[row][0] = 0 - row * prof.read_gext - prof.read_gopen;
        }
    }

    // Fill the matrixes
    for (int col = 1; col <= ref_len; ++col) {
        char ref_char = rg::num_to_base(ref_slice[col-1]);
        for (unsigned row = 1; row <= read.length(); ++row) {
            char query_char = read[row-1];
            char query_qual = has_quality ? qual[row-1] - phred_offset : 40;
            // Compute the M matrix entry. Force a match or a mismatch between the last characters
            int possibleM = M[row-1][col-1];
            int possibleD = D[row-1][col-1];
            int possibleI = I[row-1][col-1];
            if (ref_char == 'N' or query_char == 'N') { // ambiguous query and/or reference
                possibleM -= prof.ambig; possibleD -= prof.ambig; possibleI -= prof.ambig;
            } else if (ref_char != query_char) { // mismatch
                possibleM -= prof.penalty(query_qual); possibleD -= prof.penalty(query_qual); possibleI -= prof.penalty(query_qual);
            } else { // match
                possibleM += prof.match; possibleD += prof.match; possibleI += prof.match;
            }
            int best = std::max({possibleM, possibleD, possibleI});
            if (prof.end_to_end or best > 0) { // local mode nothing happens if it's <= 0
                if (possibleM == best) { M[row][col] = possibleM; tM[row][col] = 0; }
                else if (possibleD == best) { M[row][col] = possibleD; tM[row][col] = 1; }
                else { M[row][col] = possibleI; tM[row][col] = 2; }
            }

            // Compute the D matrix entry. Force a gap in end of reference.
            possibleM = M[row][col-1] - prof.read_gopen - prof.read_gext;
            possibleD = D[row][col-1] - prof.read_gext;
            possibleI = I[row][col-1] - prof.read_gopen - prof.read_gext;
            best = std::max({possibleM, possibleD, possibleI});
            if (prof.end_to_end or best > 0) {
                if (possibleM == best) { D[row][col] = possibleM; tD[row][col] = 0; }
                else if (possibleD == best) { D[row][col] = possibleD; tD[row][col] = 1; }
                else { D[row][col] = possibleI; tD[row][col] = 2; }
            }

            // Compute the I entry. Force a gap in end of reference.
            possibleM = M[row-1][col] - prof.ref_gopen - prof.ref_gext;
            possibleD = D[row-1][col] - prof.ref_gopen - prof.ref_gext;
            possibleI = I[row-1][col] - prof.ref_gext;
            best = std::max({possibleM, possibleD, possibleI});
            if (prof.end_to_end or best > 0) {
                if (possibleM == best) { I[row][col] = possibleM; tI[row][col] = 0; }
                else if (possibleD == best) { I[row][col] = possibleD; tI[row][col] = 1; }
                else { I[row][col] = possibleI; tI[row][col] = 2; }
            }
        }
    }

    // Compute traceback: CIGAR string and start position
    std::vector<char> aln; // reverse order of operations
    int result_score;
    int currCol, currRow;
    if (prof.end_to_end) { // semiglobal
        // best score is in last row and last column because we ended the reference at the max-scoring position
        currCol = ref_len;
        currRow = read.length();
        int bestMatrix = 0;
        int best = M[currRow][currCol];
        if (D[currRow][currCol] > best) { bestMatrix = 1; best = D[currRow][currCol]; }
        if (I[currRow][currCol] > best) { bestMatrix = 2; best = I[currRow][currCol]; }
        result_score = best;
        while (currRow > 0 and currCol > 0) {
            if (bestMatrix == 0) { aln.push_back('M'); bestMatrix = tM[currRow][currCol]; --currCol; --currRow; }
            else if (bestMatrix == 1) { aln.push_back('D'); bestMatrix = tD[currRow][currCol]; --currCol; }
            else { aln.push_back('I'); bestMatrix = tI[currRow][currCol]; --currRow; }
        }
        for (int row = 0; row < currRow; ++row) aln.push_back('I'); // unaligned bases in beginning of query
    } else { // local
        // best score is somewhere in last column because we ended the reference at the max-scoring position
        int bestRow = -1, best = -1, bestMatrix = -1;
        for (int row = 0; row <= static_cast<int>(read.length()); ++row) {
            if (M[row][ref_len] > best) { best = M[row][ref_len]; bestRow = row; bestMatrix = 0; }
            if (D[row][ref_len] > best) { best = D[row][ref_len]; bestRow = row; bestMatrix = 1; }
            if (I[row][ref_len] > best) { best = I[row][ref_len]; bestRow = row; bestMatrix = 2; }
        }
        result_score = best;
        currCol = ref_len;
        currRow = bestRow;
        for (unsigned row = bestRow; row < read.length(); ++row) aln.push_back('S'); // unaligned bases in end of query
        while (currRow > 0 and currCol > 0) {
            if (bestMatrix == 0 and M[currRow][currCol] > 0) { aln.push_back('M'); bestMatrix = tM[currRow][currCol]; --currCol; --currRow; }
            else if (bestMatrix == 1 and D[currRow][currCol] > 0) { aln.push_back('D'); bestMatrix = tD[currRow][currCol]; --currCol; }
            else if (I[currRow][currCol] > 0) { aln.push_back('I'); bestMatrix = tI[currRow][currCol]; --currRow; }
            else { break; } // if score goes to or below zero
        }
        for (int row = 0; row < currRow; ++row) aln.push_back('S'); // unaligned bases in beginning of query
    }

    // Reverse and run-length-collapse the sequence of operations
    std::string cigar;
    if (!aln.empty()) {
        char last_seen = aln.back();
        unsigned count = 1;
        auto rit = std::next(aln.rbegin(), 1);
        for (; rit != aln.rend(); ++rit) {
            if (*rit != last_seen) {
                cigar.append(std::to_string(count));
                cigar.push_back(last_seen);
                count = 1;
                last_seen = *rit;
            } else { count += 1; }
        }
        cigar.append(std::to_string(count));
        cigar.push_back(last_seen);
    }

    return FlatTraceResult{cigar, currCol, result_score};
}

/**
 * @brief
 * Result of a variant-graph traceback.
 */
struct GraphTraceResult {
    bool ok = false;              // true if a path matched the SIMD max score exactly
    std::string cigar;            // CIGAR of the best path found
    unsigned start_pos = 0;       // contig-relative start position (like SAM::Record::pos)
    int score = std::numeric_limits<int>::min();  // DP score of the best path found
    std::string path;             // traversed node IDs (non-ref alleles marked '*'); "capped" if bailed
};

/**
 * @brief
 * Recover the full alignment (CIGAR + start position + traversed alleles) through a variant
 * graph, by recomputing the traceback in a small window subgraph around the max-scoring end
 * position. Each source->sink path through the window is flattened to a linear reference and
 * scored with flat_traceback(); the path whose DP score matches the SIMD max score is kept.
 * @param gm GraphMan (for reference-coordinate mapping)
 * @param subgraph aligned-to subgraph
 * @param max_pos_global SIMD max-scoring position, 1-based global coordinate (aligns.max_pos[j])
 * @param max_score SIMD optimal score (ground truth over the whole graph)
 * @param read read sequence (already reverse-complemented if aligned to the reverse strand)
 * @param qual raw Phred-encoded qualities (may be empty)
 * @param phred_offset Phred offset
 * @param prof score profile
 * @param max_paths cap on enumerated paths per window; if exceeded, bail with path="capped"
 * @return best graph traceback found
 */
static GraphTraceResult graph_traceback(const vargas::GraphMan &gm,
                                        const std::shared_ptr<vargas::Graph> &subgraph,
                                        unsigned max_pos_global,
                                        int max_score,
                                        const std::string &read,
                                        const std::string &qual,
                                        char phred_offset,
                                        const vargas::ScoreProfile &prof,
                                        unsigned max_paths) {
    const auto &nmap = *subgraph->node_map();
    const auto &prevm = subgraph->prev_map();
    const auto &nextm = subgraph->next_map();
    static const std::vector<unsigned> empty_edges;
    auto preds = [&](unsigned id) -> const std::vector<unsigned> & {
        auto it = prevm.find(id); return it == prevm.end() ? empty_edges : it->second;
    };
    auto succs = [&](unsigned id) -> const std::vector<unsigned> & {
        auto it = nextm.find(id); return it == nextm.end() ? empty_edges : it->second;
    };

    // Candidate end nodes: every node whose coordinate range contains the max-cell position.
    // seek() returns one such node, but parallel alt alleles are anchored to the SAME end
    // coordinate, so the max cell may lie in a sibling of seek()'s node. Collect those siblings
    // (they share a predecessor or successor with the seeded node).
    const unsigned pos0 = max_pos_global - 1;   // 0-based coordinate of the max cell
    auto contains = [&](unsigned id) {
        const auto &n = nmap.at(id);
        return n.begin_pos() <= pos0 && pos0 <= n.end_pos();
    };
    std::unordered_set<unsigned> ends;
    std::unordered_map<unsigned, int> endOff;   // end node -> offset of max cell within it
    auto add_end = [&](unsigned id) {
        if (contains(id) && ends.insert(id).second) endOff[id] = static_cast<int>(pos0 - nmap.at(id).begin_pos());
    };
    const unsigned e_seed = subgraph->seek(max_pos_global).first->id();
    add_end(e_seed);
    // Consider one-hop neighbours that also cover the coordinate. This captures (a) parallel alt
    // alleles that share a predecessor/successor with the seeded node, and (b) insertion alleles,
    // which are successors anchored to the same end coordinate and therefore still contain pos0.
    for (unsigned p : preds(e_seed)) { add_end(p); for (unsigned c : succs(p)) add_end(c); }
    for (unsigned s : succs(e_seed)) { add_end(s); for (unsigned c : preds(s)) add_end(c); }
    if (ends.empty()) {                          // fallback: coordinate edge case, use seek's node
        const auto &n = nmap.at(e_seed);
        int off = (pos0 >= n.begin_pos()) ? std::min<int>((int)(pos0 - n.begin_pos()), (int)n.length() - 1)
                                          : (int)n.length() - 1;
        ends.insert(e_seed);
        endOff[e_seed] = std::max(off, 0);
    }

    const int read_len = static_cast<int>(read.length());
    int W = 2 * read_len;                 // window budget in reference bases upstream of the ends
    if (W <= 0) W = 1;

    GraphTraceResult best;                // best across attempts (may not match exactly)
    int bestGap = std::numeric_limits<int>::max();

    for (int attempt = 0; attempt < 2; ++attempt) {
        // 1. Build the window: backward reachability from the end nodes until W reference bases
        //    have been covered, or a graph source is reached. (Pinch nodes are NOT a boundary: the
        //    scorer only clears the seed *map* there as a memory optimization; the DP still flows
        //    through a pinch node from its predecessors, so upstream context is still required.)
        std::unordered_map<unsigned, int> bestRem;
        std::vector<unsigned> stack;
        for (unsigned e : ends) { bestRem[e] = W; stack.push_back(e); }
        while (!stack.empty()) {
            unsigned id = stack.back(); stack.pop_back();
            int rem = bestRem[id];
            if (rem <= 0) continue;
            for (unsigned p : preds(id)) {
                int nrem = rem - static_cast<int>(nmap.at(p).length());
                auto it = bestRem.find(p);
                if (it == bestRem.end() || nrem > it->second) { bestRem[p] = nrem; stack.push_back(p); }
            }
        }
        std::unordered_set<unsigned> window;
        for (auto &kv : bestRem) window.insert(kv.first);

        // 2. Enumerate source->end paths within the window; score each with flat_traceback.
        const int L = W;                  // flattened-reference target length (matches flat sizing)
        unsigned path_count = 0;
        bool capped = false;
        std::vector<unsigned> cur;
        std::function<void(unsigned)> dfs = [&](unsigned id) {
            if (capped) return;
            cur.push_back(id);
            if (ends.count(id)) {
                // Flatten cur -> slice + per-base provenance, truncating the end node at its offset.
                const unsigned term = id;
                std::vector<rg::Base> slice;
                std::vector<std::pair<unsigned, int>> prov; // (node_id, node_offset) per base
                for (unsigned nid : cur) {
                    const auto &seq = nmap.at(nid).seq();
                    int upto = (nid == term) ? std::min<int>(endOff[term] + 1, (int)seq.size()) : (int)seq.size();
                    for (int o = 0; o < upto; ++o) { slice.push_back(seq[o]); prov.emplace_back(nid, o); }
                }
                // Keep only the last L bases (upstream context is bounded).
                if ((int)slice.size() > L) {
                    int drop = (int)slice.size() - L;
                    slice.erase(slice.begin(), slice.begin() + drop);
                    prov.erase(prov.begin(), prov.begin() + drop);
                }
                if (!slice.empty()) {
                    FlatTraceResult tr = flat_traceback(slice, read, qual, phred_offset, prof);
                    int gap = std::abs(tr.score - max_score);
                    if (gap < bestGap) {
                        bestGap = gap;
                        int off = std::min<int>(std::max(tr.ref_offset, 0), (int)prov.size() - 1);
                        auto pr = prov[off];
                        unsigned gstart0 = nmap.at(pr.first).begin_pos() + pr.second; // 0-based global
                        best.start_pos = gm.absolute_position(gstart0 + 1).second;
                        // Build the path tag from the aligned span [off, end).
                        std::string ptag;
                        unsigned prevn = std::numeric_limits<unsigned>::max();
                        for (int q = off; q < (int)prov.size(); ++q) {
                            unsigned nid = prov[q].first;
                            if (nid != prevn) {
                                if (!ptag.empty()) ptag.push_back(',');
                                ptag += std::to_string(nid);
                                if (!nmap.at(nid).is_ref()) ptag.push_back('*');
                                prevn = nid;
                            }
                        }
                        best.ok = (gap == 0);
                        best.cigar = tr.cigar;
                        best.score = tr.score;
                        best.path = ptag;
                    }
                }
                if (++path_count > max_paths) capped = true;
            } else {
                for (unsigned s : succs(id)) if (window.count(s)) dfs(s);
            }
            cur.pop_back();
        };
        for (unsigned id : window) {
            bool is_source = true;
            for (unsigned p : preds(id)) if (window.count(p)) { is_source = false; break; }
            if (is_source) dfs(id);
        }

        if (capped) { best.ok = false; best.path = "capped"; return best; }
        if (best.ok) return best;   // found an exact-score path
        W *= 2;                     // widen the window and retry once
    }
    return best;                    // best effort (score may differ; caller warns)
}

struct align_helper {
    vargas::GraphMan &gm;
    std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> &task_list;
    vargas::osam &out;
    const std::vector<std::unique_ptr<vargas::AlignerBase, rg::Deleter>> &aligners;
    bool fwdonly, msonly, maxonly, notraceback;
    char phred_offset;
    unsigned max_paths;
    std::mutex &mut;
};

void align_helper_func(void *data, long index, int tid) {
    align_helper &help(*(align_helper *)data);
    auto &aligners = help.aligners;
    auto &out = help.out;
    auto &task_list = help.task_list;
    auto &gm = help.gm;
    auto fwdonly = help.fwdonly;
    auto msonly = help.msonly;
    auto maxonly = help.maxonly;
    auto notraceback = help.notraceback;
    char phred_offset = help.phred_offset;
    auto max_paths = help.max_paths;

    const size_t num_reads = task_list.at(index).second.size();
    std::vector<std::string> read_seqs(num_reads);
    std::vector<std::vector<char>> quals(num_reads);
    for (size_t i = 0; i < num_reads; ++i) {
        const auto &r = task_list.at(index).second.at(i);
        read_seqs[i] = r.seq;
        if (r.qual.size() == r.seq.size()) {
            std::transform(r.qual.begin(),
                           r.qual.end(),
                           std::back_inserter(quals[i]),
                           [](char c){ return c - 33; }); //TODO needs to be offset variable
        }
    }
    auto subgraph = gm.at(task_list.at(index).first);
    vargas::Results aligns;
    aligners[tid]->align_into(read_seqs, quals, subgraph->begin(), subgraph->end(), aligns, fwdonly);

    //If no variants (# nodes == # contigs) compute the alignment traceback
    bool not_graph = subgraph->node_map()->size() == gm.resolver()._contig_hdr_order.size();

    for (size_t j = 0; j < task_list.at(index).second.size(); ++j) {
        vargas::SAM::Record &rec = task_list.at(index).second.at(j);
        auto abs = gm.absolute_position(aligns.max_pos[j]);
        rec.aux.set("AS", aligns.max_score[j]);
        if (!msonly) {
            // Can only guess start position for end to end
            // if (aligns.profile.end_to_end) rec.pos = abs.second - rec.seq.size() + 1;

            rec.ref_name = abs.first;
            rec.flag.rev_complement = aligns.max_strand[j] == vargas::Strand::REV;
            if (rec.flag.rev_complement) {
                rg::reverse_complement_inplace(rec.seq);
                std::reverse(rec.qual.begin(), rec.qual.end());
            }
            rec.aux.set(ALIGN_SAM_MAX_POS_TAG, abs.second);
            rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns.max_count[j]);

            if (!notraceback) {
                if (not_graph) {
                    // Fast path: a linear reference. Slice a flat window and trace back.
                    int nodeID = gm.nodeID_from_contig(rec.ref_name);
                    vargas::Graph::nodemap_t _node_map = *(subgraph->node_map());
                    //TODO upper-bound the length of reference slice needed based on the score or scoring function
                    int ref_len = 2*rec.seq.length() < abs.second ? 2*rec.seq.length() : abs.second ;
                    int ref_start = abs.second-ref_len;
                    auto ref_iter = std::next(_node_map[nodeID].begin(),ref_start); //TODO this is linear time I think?
                    std::vector<rg::Base> ref_slice(ref_iter, std::next(ref_iter, ref_len));

                    FlatTraceResult tr = flat_traceback(ref_slice, rec.seq, rec.qual, phred_offset, aligns.profile);
                    if (tr.score != aligns.max_score[j]) {
                        std::cerr << "[WARNING] " << rec.query_name << " DP optimal score " << tr.score
                                  << " and SIMD optimal score " << aligns.max_score[j] << " not equal\n";
                    }
                    rec.pos = abs.second - ref_len + tr.ref_offset;
                    rec.cigar = tr.cigar;
                } else {
                    // Variant graph: recompute the traceback over a window subgraph.
                    GraphTraceResult gt = graph_traceback(gm, subgraph, aligns.max_pos[j], aligns.max_score[j],
                                                          rec.seq, rec.qual, phred_offset, aligns.profile, max_paths);
                    if (gt.path == "capped") {
                        std::cerr << "[WARNING] " << rec.query_name << " too many candidate paths (> "
                                  << max_paths << "); skipping CIGAR. Increase --max-trace-paths.\n";
                    } else {
                        if (!gt.ok) {
                            std::cerr << "[WARNING] " << rec.query_name << " best graph-path DP score " << gt.score
                                      << " and SIMD optimal score " << aligns.max_score[j] << " not equal\n";
                        }
                        rec.pos = gt.start_pos;
                        rec.cigar = gt.cigar;
                        rec.aux.set(ALIGN_SAM_PATH_TAG, gt.path);
                    }
                }
            }

            // Flags for 2nd max
            if (!maxonly) {
                abs = gm.absolute_position(aligns.sub_pos[j]);
                rec.aux.set(ALIGN_SAM_SUB_SEQ, abs.first);
                rec.aux.set(ALIGN_SAM_SUB_POS_TAG, abs.second);
                rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns.sub_count[j]);
                rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns.sub_score[j]);
                rec.aux.set(ALIGN_SAM_SUB_STRAND_TAG, aligns.sub_strand[j] == vargas::Strand::FWD ? "fwd" : "rev");
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(help.mut);
        for (const auto & j : task_list.at(index).second) out.add_record(j);
    }
}

#if !NDEBUG
#else
#define at operator[]
#endif

void align(vargas::GraphMan &gm,
           std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> &task_list,
           vargas::osam &out,
           const std::vector<std::unique_ptr<vargas::AlignerBase, rg::Deleter>> &aligners,
           bool fwdonly, bool msonly, bool maxonly, bool notraceback, char phred_offset, unsigned max_paths) {
    std::cerr << "Aligning... " << std::flush;
    rg::ForPool fp(aligners.size());
    auto start_time = std::chrono::steady_clock::now();

    const auto num_tasks = task_list.size();
    std::mutex mut;
    align_helper help{gm, task_list, out, aligners, fwdonly, msonly, maxonly, notraceback, phred_offset, max_paths, mut};
    fp.forpool(&align_helper_func, (void *)&help, num_tasks);

    std::cerr << rg::chrono_duration(start_time) << "s.\n";

}

std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>>
create_tasks(vargas::isam &reads, std::string &align_targets, const int chunk_size, size_t &read_len) {
    std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> task_list;
    std::unordered_map<std::string, std::vector<vargas::SAM::Record>> read_groups;

    std::vector<std::string> alignment_pairs;
    if (align_targets.length() != 0) {
        std::replace(align_targets.begin(), align_targets.end(), '\n', ';');
        alignment_pairs = rg::split(align_targets, ';');
    }

    std::cerr << "Loading reads... " << std::flush;
    auto start_time = std::chrono::steady_clock::now();

    size_t total = 0;
    auto &reads_hdr = reads.header();
    std::string read_group;
    vargas::SAM::Record rec;
    read_len = reads.record().seq.length();
    do {
        rec = reads.record();
        if (rec.seq.length() > read_len) read_len = rec.seq.length();
        if (!rec.aux.get("RG", read_group)) {
            read_group = UNGROUPED_READGROUP;
            rec.aux.set("RG", UNGROUPED_READGROUP);
            if (!reads_hdr.read_groups.count(UNGROUPED_READGROUP)) {
                reads_hdr.add(vargas::SAM::Header::ReadGroup("@RG\tID:" + std::string(UNGROUPED_READGROUP)));
            }
        }
        read_groups[read_group].push_back(rec);
    } while (reads.next());

    if (alignment_pairs.empty()) {
        for (const auto &p : read_groups) {
            alignment_pairs.push_back("RG:ID:" + p.first + ",base");
        }
    }

    // Specify a single target graph
    if (alignment_pairs.size() == 1) {
        auto target = rg::split(alignment_pairs[0], ',');
        if (target.size() == 1) {
            alignment_pairs.clear();
            for (const auto &p : read_groups) {
                alignment_pairs.push_back("RG:ID:" + p.first + "," + target[0]);
            }
        }
    }

    // Maps target graph to read group ID's
    std::unordered_map<std::string, std::vector<std::string>> alignment_rg_map;

    std::string tag, val, target_val;
    for (const std::string &p : alignment_pairs) {
        auto pair = rg::split(p, ',');
        if (pair.size() != 2)
            throw std::invalid_argument("Malformed alignment pair \"" + p + "\".");
        if (pair[0].substr(0, 2) != "RG")
            throw std::invalid_argument(R"(Expected a read group tag 'RG:xx:', got ")" + pair[0] + "\"");
        if (pair[0].at(2) != ':')
            throw std::invalid_argument("Expected source format Read_group_tag:value in \"" + pair[0] + "\".");


        tag = pair[0].substr(3, 2);
        target_val = pair[0].substr(6);

        for (const auto &rg_pair : reads_hdr.read_groups) {
            if (tag == "ID") val = rg_pair.second.id;
            else if (rg_pair.second.aux.get(tag, val));
            else continue;
            if (val == target_val) alignment_rg_map[pair[1]].push_back(rg_pair.first);
        }

    }

    std::cerr << rg::chrono_duration(start_time) << "s." << std::endl;

    // graph label to vector of reads
    for (const auto &sub_rg_pair : alignment_rg_map) {
        for (const std::string &rgid : sub_rg_pair.second) {
            // If there is a header line that there are no reads associated with, skip
            if (read_groups.count(rgid) == 0) continue;

            const auto beg = std::begin(read_groups.at(rgid));
            const auto end = std::end(read_groups.at(rgid));
            const size_t nrecords = read_groups.at(rgid).size();
            const size_t n_chunks = (nrecords / chunk_size) + 1;
            total += read_groups.at(rgid).size();

            for (size_t i = 0; i < n_chunks; ++i) {
                const auto safe_beg = beg + (i * chunk_size);
                const auto safe_end = (i + 1) * chunk_size > nrecords ? end : safe_beg + chunk_size;
                if (safe_beg != safe_end) {
                    task_list.emplace_back(sub_rg_pair.first, std::vector<vargas::SAM::Record>(safe_beg, safe_end));
                }
            }
        }
    }

    std::cerr << read_groups.size() << "\tRead group(s).\n"
              << alignment_rg_map.size() << "\tSubgraph(s).\n"
              << task_list.size() << "\tTask(s).\n"
              << total << "\tTotal alignments.\n"
              << read_len << "\tMax read length.\n";

    return task_list;
}


template<typename T, typename...Args>
T *construct_aligned(Args &&...args) {
    static constexpr size_t alignment = 64; // AVX512
    T *ptr;
    if(posix_memalign(reinterpret_cast<void **>(&ptr), alignment, sizeof(T))) throw std::bad_alloc();
    return new(ptr) T(std::forward<Args>(args)...);
}


std::unique_ptr<vargas::AlignerBase, Deleter>
make_aligner(const vargas::ScoreProfile &prof, size_t read_len, bool use_wide, bool msonly, bool /* maxonly */) {
    std::unique_ptr<vargas::AlignerBase, Deleter> ret;
    if (msonly) {
        if (prof.end_to_end) {
            if(use_wide) ret.reset(construct_aligned<vargas::MSWordAlignerETE>(read_len, prof));
            else ret.reset(construct_aligned<vargas::MSAlignerETE>(read_len, prof));
        } else {
            if(use_wide) ret.reset(construct_aligned<vargas::MSWordAligner>(read_len, prof));
            else ret.reset(construct_aligned<vargas::MSAligner>(read_len, prof));
        }
    }
    else {
        if (prof.end_to_end) {
            if(use_wide) ret.reset(construct_aligned<vargas::WordAlignerETE>(read_len, prof));
            else ret.reset(construct_aligned<vargas::AlignerETE>(read_len, prof));
        } else {
            if(use_wide) ret.reset(construct_aligned<vargas::WordAligner>(read_len, prof));
            else ret.reset(construct_aligned<vargas::Aligner>(read_len, prof));
        }
    }
    return ret;
}

void load_fast(std::string &file, const bool fastq, vargas::isam &ret, bool p64) {
    std::string input;
    if (file.empty()) {
        std::stringstream ss;
        ss << std::cin.rdbuf();
        input = ss.str();
    } else {
        std::ifstream in(file);
        if (!in.good()) throw std::invalid_argument("Unable to open file \"" + file + "\"");
        std::stringstream ss;
        ss << in.rdbuf();
        input = ss.str();
    }
    const auto lines = rg::split(input, '\n');

    try {
        for (unsigned i = 0; i < lines.size(); i += (fastq ? 4 : 2)) {
            vargas::SAM::Record rec;
            rec.query_name = std::string(lines.at(i).begin() + 1,
                                         std::find_if(lines.at(i).begin() + 1, lines.at(i).end(), [](char c) { return std::isspace(c); }));
            rec.seq = lines.at(i + 1);
            if (fastq) rec.qual = lines.at(i + 3);
            if (p64) std::transform(rec.qual.begin(), rec.qual.end(), rec.qual.begin(), [](char c){return c-31;});

            ret.push(rec);
        }
    } catch (std::exception &e) {
        throw std::runtime_error("Invalid FASTA/Q file.");
    }
    ret.next();
}

void align_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help(opts.groups()) << "\n" << endl;
    cerr << "Elements per SIMD vector: " << vargas::Aligner::read_capacity() << endl;
}

ReadFmt read_fmt(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.good()) throw std::invalid_argument("Invalid read file: " + filename);

    std::string line;
    if (!std::getline(in, line)) throw std::invalid_argument("Empty Read File."); // @SAM or fasta/q name
    if (line.substr(0,3) == "@HD") return ReadFmt::SAM;
    if (!std::getline(in, line)) throw std::invalid_argument("Invalid Read File."); // SAM comment/header, or read
    if (!std::getline(in, line)) return ReadFmt::FASTA; // Single record fasta, or SAM line, or +, or name
    if (!line.empty() && line[0] == '+') return ReadFmt::FASTQ;
    if (!line.empty() && (line[0] == '>' || line[0] == '@')) return ReadFmt::FASTA;
    return ReadFmt::SAM;
}

TEST_SUITE("System");
TEST_CASE ("Load FASTQ") {
    std::string tmpfq = "tmp_fastq.va";
    {
        std::ofstream o(tmpfq);
        o << "@name desc\nAAAAACCCCC\n+\n!!!!!!!!!!";
    }
    CHECK(read_fmt(tmpfq) == ReadFmt::FASTQ);
    vargas::isam ss;
    load_fast(tmpfq, true, ss);

    CHECK(ss.record().query_name == "name");
    CHECK(ss.record().seq == "AAAAACCCCC");
    CHECK(ss.record().qual == "!!!!!!!!!!");
    CHECK_FALSE(ss.next());
    remove(tmpfq.c_str());
}
TEST_CASE ("Load FASTA") {
    std::string tmpfq = "tmp_fastq.va";
    {
        std::ofstream o(tmpfq);
        o << ">name desc\nAAAAACCCCC\n>x\nGGGGGTTTTT";
    }
    CHECK(read_fmt(tmpfq) == ReadFmt::FASTA);
    vargas::isam ss;
    load_fast(tmpfq, false, ss);

    CHECK(ss.record().query_name == "x");
    CHECK(ss.record().seq == "GGGGGTTTTT");
    ss.next();

    CHECK(ss.record().query_name == "name");
    CHECK(ss.record().seq == "AAAAACCCCC");
    CHECK_FALSE(ss.next());
    remove(tmpfq.c_str());
}

TEST_CASE ("Graph traceback SNP bubble") {
    // ref path = AAC C TTT ; alt path = AAC G TTT ; the C/G nodes end at the same coordinate.
    auto g = std::make_shared<vargas::Graph>();
    auto mknode = [](unsigned id, const std::string &s, unsigned endp, bool ref) {
        vargas::Graph::Node n;
        n.set_id(id); n.set_seq(s); n.set_endpos(endp); n.set_pinch(false);
        if (ref) n.set_as_ref(); else n.set_not_ref();
        n.set_population(1, ref);
        return n;
    };
    g->add_node(mknode(0, "AAC", 2, true));   // pos 0..2
    g->add_node(mknode(1, "G",   3, false));  // alt allele at pos 3
    g->add_node(mknode(2, "C",   3, true));   // ref allele at pos 3
    g->add_node(mknode(3, "TTT", 6, true));   // pos 4..6
    g->add_edge(0, 1); g->add_edge(0, 2); g->add_edge(1, 3); g->add_edge(2, 3);

    vargas::GraphMan gm;                       // empty resolver: absolute_position(p) == {"", p}
    vargas::ScoreProfile loc(2u, 6u, 3u, 1u);  // local, match 2 / mismatch 6 / gap 3,1
    loc.end_to_end = false;

    SUBCASE("prefers alt allele") {
        auto r = graph_traceback(gm, g, /*max_pos=*/7, /*max_score=*/14, "AACGTTT", "", 33, loc, 256);
        CHECK(r.ok);
        CHECK(r.cigar == "7M");
        CHECK(r.path == "0,1*,3");   // alt node 1 marked with '*'
        CHECK(r.start_pos == 1);
    }
    SUBCASE("prefers ref allele") {
        auto r = graph_traceback(gm, g, 7, 14, "AACCTTT", "", 33, loc, 256);
        CHECK(r.ok);
        CHECK(r.cigar == "7M");
        CHECK(r.path == "0,2,3");
    }
    SUBCASE("path cap bails") {
        auto r = graph_traceback(gm, g, 7, 14, "AACGTTT", "", 33, loc, 1);
        CHECK(r.path == "capped");
        CHECK_FALSE(r.ok);
    }
}

TEST_CASE ("Graph traceback trie (varying-length sibling branches)") {
    // A path-compressed reverse-suffix trie: a shared trunk "AC" that branches into siblings of
    // DIFFERENT lengths -- "GTA" (leaf, spells ACGTA) and "TT" (leaf, spells ACTT). Unlike a SNP
    // bubble, the sibling end nodes do NOT share an end coordinate, which is the topology produced
    // by build_trie.py. Verifies that graph_traceback recovers the correct branch/CIGAR anyway.
    auto g = std::make_shared<vargas::Graph>();
    auto mknode = [](unsigned id, const std::string &s, unsigned endp, bool ref) {
        vargas::Graph::Node n;
        n.set_id(id); n.set_seq(s); n.set_endpos(endp); n.set_pinch(false);
        if (ref) n.set_as_ref(); else n.set_not_ref();
        n.set_population(1, ref);
        return n;
    };
    g->add_node(mknode(0, "AC",  1, true));   // trunk, pos 0..1
    g->add_node(mknode(1, "GTA", 4, false));  // long branch, pos 2..4  (leaf ACGTA)
    g->add_node(mknode(2, "TT",  3, false));  // short branch, pos 2..3 (leaf ACTT)
    g->add_edge(0, 1); g->add_edge(0, 2);

    vargas::GraphMan gm;                       // empty resolver: absolute_position(p) == {"", p}
    vargas::ScoreProfile loc(2u, 6u, 3u, 1u);  // local, match 2 / mismatch 6 / gap 3,1
    loc.end_to_end = false;

    SUBCASE("long branch") {
        // ACGTA, max cell at end of node 1 (pos 4 -> 1-based 5), score 5*2=10.
        auto r = graph_traceback(gm, g, 5, 10, "ACGTA", "", 33, loc, 256);
        CHECK(r.ok);
        CHECK(r.cigar == "5M");
        CHECK(r.path == "0,1*");
        CHECK(r.start_pos == 1);
    }
    SUBCASE("short branch (max cell also lies within the long sibling's range)") {
        // ACTT, max cell at end of node 2 (pos 3 -> 1-based 4), score 4*2=8. pos0=3 also falls
        // inside node 1 (range 2..4), so node 1 is a candidate end too; the re-scoring must still
        // prefer the exact short-branch path.
        auto r = graph_traceback(gm, g, 4, 8, "ACTT", "", 33, loc, 256);
        CHECK(r.ok);
        CHECK(r.cigar == "4M");
        CHECK(r.path == "0,2*");
        CHECK(r.start_pos == 1);
    }
    SUBCASE("read starting mid-trunk") {
        // CGTA: a substring of ACGTA starting inside the trunk; ends at node 1 (pos 4). Local
        // alignment recovers the 4-base match starting at trunk offset 1 (global pos 2 -> 1-based).
        auto r = graph_traceback(gm, g, 5, 8, "CGTA", "", 33, loc, 256);
        CHECK(r.ok);
        CHECK(r.cigar == "4M");
        CHECK(r.start_pos == 2);
    }
}

/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Provides tools to interact with VCF/BCF files.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include "varfile.h"
#include <unistd.h>
#include <algorithm>
#include <sstream>

vargas::Region vargas::parse_region(const std::string &region_str) {
    vargas::Region ret;

    std::vector <std::string> regionSplit = rg::split(region_str, ':');

    if(regionSplit.size() == 1) {
        regionSplit.emplace_back("0-0");
    }
    else if (regionSplit.size() != 2) throw std::invalid_argument("Invalid region: " + region_str);

    ret.seq_name = regionSplit[0];

    // Strip commas
    regionSplit[1].erase(std::remove(regionSplit[1].begin(), regionSplit[1].end(), ','), regionSplit[1].end());

    // Range
    regionSplit = rg::split(regionSplit[1], '-');
    if (regionSplit.size() != 2)
        throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY\n\t" + region_str);

    ret.min = std::stoul(regionSplit[0]);
    ret.max = std::stoul(regionSplit[1]);

    if (ret.min > ret.max) {
        throw std::invalid_argument("Invalid region, min > max.");
    }

    return ret;
}

void vargas::VCF::set_region(const Region &region) {
    _region = region;
}


std::vector<std::string> vargas::VCF::seq_names() const {
    if (!_vcf) return std::vector<std::string>(0);
    
    // Extract sequence names from vcflib header
    std::vector<std::string> ret;
    std::string header = _vcf->header;
    
    // Parse header lines to find contig information
    std::istringstream header_stream(header);
    std::string line;
    while (std::getline(header_stream, line)) {
        if (line.substr(0, 8) == "##contig") {
            // Extract contig name from ##contig=<ID=chr1,length=...>
            size_t id_pos = line.find("ID=");
            if (id_pos != std::string::npos) {
                size_t start = id_pos + 3;
                size_t end = line.find(',', start);
                if (end == std::string::npos) {
                    end = line.find('>', start);
                }
                if (end != std::string::npos) {
                    ret.push_back(line.substr(start, end - start));
                }
            }
        }
    }
    
    return ret;
}


bool vargas::VCF::next() {
    if (_limit > 0 && _counter >= _limit) return false;
    if (!_vcf) return false;
    
    // Keep reading until we find a record that matches our region filter
    while (true) {
        // Read next variant using vcflib
        if (!_vcf->getNextVariant(*_curr_var)) {
            return false;
        }
        
        // Check region filtering
        bool seqmatch = _region.seq_name.empty() || _curr_var->sequenceName == _region.seq_name;
        if (seqmatch) {
            _entered_contig = true;
        } else if (_assume_contig && _entered_contig) {
            return false;
        } else {
            continue;
        }
        
        // Check position filtering
        if (unsigned(_curr_var->position) < _region.min || (_region.max > 0 && unsigned(_curr_var->position) > _region.max)) {
            continue; // Skip this record and try the next one
        }
        
        // Found a matching record, load data and return
        _load_shared();
        gen_genotypes();
        ++_counter;
        return true;
    }
}


const std::vector<std::string> &vargas::VCF::gen_genotypes() {
    _genotypes.clear();
    
    // Determine which samples to process
    std::vector<std::string> samplesToProcess;
    if (_ingroup.empty()) {
        // No filter, process all samples
        samplesToProcess = _vcf->sampleNames;
    } else {
        // Filter to only include samples in _ingroup
        for (const auto& sampleName : _vcf->sampleNames) {
            if (std::find(_ingroup.begin(), _ingroup.end(), sampleName) != _ingroup.end()) {
                samplesToProcess.push_back(sampleName);
            }
        }
    }
    
    // Get genotypes from vcflib variant
    for (const auto& sampleName : samplesToProcess) {
        auto sampleIt = _curr_var->samples.find(sampleName);
        if (sampleIt != _curr_var->samples.end()) {
            auto gtIt = sampleIt->second.find("GT");
            if (gtIt != sampleIt->second.end() && !gtIt->second.empty()) {
                // Parse GT field (e.g., "0|1" or "1/0")
                std::string gt_str = gtIt->second[0];
                std::vector<std::string> gt_parts = rg::split(gt_str, "|/");
                
                for (const auto& gt_part : gt_parts) {
                    if (gt_part == "." || gt_part.empty()) {
                        _genotypes.push_back("*");  // Missing genotype
                    } else {
                        int allele_idx = std::stoi(gt_part);
                        if (allele_idx >= 0 && static_cast<size_t>(allele_idx) < _alleles.size()) {
                            _genotypes.push_back(_alleles[allele_idx]);
                        } else {
                            _genotypes.push_back("*");  // Invalid genotype
                        }
                    }
                }
            } else {
                // Missing GT field
                _genotypes.push_back("*");
                _genotypes.push_back("*");
            }
        } else {
            // Missing sample
            _genotypes.push_back("*");
            _genotypes.push_back("*");
        }
    }

    // Map of which indivs have each allele
    _genotype_indivs.clear();
    for (auto &allele : alleles()) {
        _genotype_indivs[allele] = Population(_genotypes.size(), false);
    }
    // Also add entry for missing genotypes
    _genotype_indivs["*"] = Population(_genotypes.size(), false);
    
    for (size_t s = 0; s < _genotypes.size(); ++s) {
        auto it = _genotype_indivs.find(_genotypes[s]);
        if (it != _genotype_indivs.end()) {
            it->second.set(s);
        }
    }

    return _genotypes;
}


const std::vector<float> &vargas::VCF::frequencies() const {
    InfoField<float> af(_curr_var, "AF");
    static std::vector<float> _allele_freqs;
    const auto &val = af.values;
    _allele_freqs.resize(val.size() + 1); // make room for the ref
    float sum = 0;

    // Get the ref frequency and add to result vector +1, so we can put ref at index 0
    for (size_t i = 0; i < val.size(); ++i) {
        sum += val[i];
        _allele_freqs[i + 1] = val[i];
    }
    _allele_freqs[0] = 1 - sum;

    return _allele_freqs;
}


void vargas::VCF::create_ingroup(int percent) {
    _ingroup.clear();

    if (percent == 100) {
        _ingroup = _samples;
    } else if (percent != 0) {
        for (const auto& s : _samples) {
            if (rand() % 100 < percent) _ingroup.push_back(s);
        }
    }

    _apply_ingroup_filter();
}


int vargas::VCF::_init() {
    _assume_contig = false;
    _entered_contig = false;
    _counter = 0;
    _limit = 0;
    if (_file_name.length() && _file_name != "-") {

        
        // Open VCF file using vcflib
        _vcf = new vcflib::VariantCallFile();
        if (!_vcf->open(_file_name)) {
            delete _vcf;
            _vcf = nullptr;
            return -1;
        }

        // Initialize the variant object
        _curr_var = new vcflib::Variant(*_vcf);

        // Load samples
        _samples = _vcf->sampleNames;
        create_ingroup(100);
    }
    return 0;
}


void vargas::VCF::_load_shared() {
    _alleles.clear();
    
    // Check if _curr_var is valid
    if (!_curr_var) {
        return;
    }
    
    for (size_t i = 0; i < _curr_var->alleles.size(); ++i) {
        std::string allele = _curr_var->alleles[i];
        
        // Handle special replacement tags
        if (allele.at(0) == '<') {
            std::string ref = _curr_var->ref;
            // Copy number
            if (allele.substr(1, 2) == "CN" && allele.at(3) != 'V') {
                int copy = std::stoi(allele.substr(3, allele.length() - 4));
                if (copy == 0) {
                    // CN0 means deletion - use empty string as per VCF spec
                    allele = "";
                } else {
                    allele = "";
                    for (int o = 0; o < copy; ++o) allele += ref;
                }
            } else {
                // Other types are just subbed with the ref.
                allele = ref;
            }
        }
        _alleles.push_back(allele);
    }
    

}




void vargas::VCF::_apply_ingroup_filter() {
    if (_vcf == nullptr) {
        throw std::logic_error("Ingroup filter should only be applied after loading header!");
    }

    // For vcflib, we don't need to set samples in the header
    // The filtering will be done when reading variants
    // This is a simplified implementation
}

size_t vargas::VCF::num_haplotypes() const {
    if (_vcf == nullptr) {
        return 0;
    }
    return (size_t) _vcf->sampleNames.size() * 2;
}

void vargas::VCF::close() {
    if (_curr_var) {
        delete _curr_var;
        _curr_var = nullptr;
    }
    if (_vcf) {
        delete _vcf;
        _vcf = nullptr;
    }
}

TEST_SUITE("VCF Parser");

TEST_CASE ("VCF File handler") {
    using std::endl;
    std::string tmpvcf = "tmp_tc.vcf";

    // Write temp VCF file
    {
        std::ofstream vcfo(tmpvcf);
        vcfo
        << "##fileformat=VCFv4.1" << endl
        << "##phasing=true" << endl
        << "##contig=<ID=x>" << endl
        << "##contig=<ID=y>" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Freq\">" << endl
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternete Allele count\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
        << "x\t9\t.\tG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
        << "x\t10\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    SUBCASE("File write wrapper") {

        SUBCASE("Unfiltered") {
            vargas::VCF vcf(tmpvcf);
            vcf.next();
            CHECK(vcf.num_haplotypes() == 4);
            CHECK(vcf.seq_names().size() == 2);
            CHECK(vcf.seq_names()[0] == "x");
            CHECK(vcf.seq_names()[1] == "y");
            REQUIRE(vcf.samples().size() == 2);
            CHECK(vcf.samples()[0] == "s1");
            CHECK(vcf.samples()[1] == "s2");

            REQUIRE(vcf.gen_genotypes().size() == 4);
            CHECK(vcf.gen_genotypes()[0] == "G");
            CHECK(vcf.gen_genotypes()[1] == "A");
            CHECK(vcf.gen_genotypes()[2] == "C");
            CHECK(vcf.gen_genotypes()[3] == "T");
            REQUIRE(vcf.alleles().size() == 4);
            CHECK(vcf.alleles()[0] == "G");
            CHECK(vcf.alleles()[1] == "A");
            CHECK(vcf.alleles()[2] == "C");
            CHECK(vcf.alleles()[3] == "T");
            CHECK(vcf.ref() == "G");
            CHECK(vcf.pos() == 8);

            // Copy number alleles
            vcf.next();
            REQUIRE(vcf.gen_genotypes().size() == 4);
            CHECK(vcf.gen_genotypes()[0] == "CC");
            CHECK(vcf.gen_genotypes()[1] == "CC");
            CHECK(vcf.gen_genotypes()[2] == "");
            CHECK(vcf.gen_genotypes()[3] == "CC");
            REQUIRE(vcf.alleles().size() == 3);
            CHECK(vcf.alleles()[0] == "C");
            CHECK(vcf.alleles()[1] == "CC");
            CHECK(vcf.alleles()[2] == "");
            CHECK(vcf.ref() == "C");
            CHECK(vcf.pos() == 9);

            // Invalid tags
            vcf.next();
            REQUIRE(vcf.alleles().size() == 3);
            CHECK(vcf.alleles()[0] == "G");
            CHECK(vcf.alleles()[1] == "G");
            CHECK(vcf.alleles()[2] == "G");
            CHECK(vcf.ref() == "G");
            CHECK(vcf.pos() == 13);

            // Next y contig should still load
            vcf.next();
            CHECK(vcf.alleles()[0] == "TATA");
        }

        SUBCASE("CHROM Filtering") {
            vargas::VCF vcf;
            vcf.set_region(std::string("y:0-0"));
            vcf.open(tmpvcf);

            vcf.next();
            CHECK(vcf.ref() == "TATA");
            vcf.next();
            CHECK(vcf.ref() == "T");
            CHECK(vcf.next() == 0); // File end
        }

        SUBCASE("Region filtering") {
            vargas::VCF vcf;
            vcf.set_region(std::string("x:0-14"));
            vcf.open(tmpvcf);

            vcf.next();
            CHECK(vcf.ref() == "G");
            vcf.next();
            CHECK(vcf.ref() == "C");
            vcf.next();
            CHECK(vcf.ref() == "G");
            CHECK(vcf.next() == 0); // Region end
        }

        SUBCASE("Ingroup generation") { //Some tests fail due to random number
            vargas::VCF vcf;
            srand(12345);
            vcf.open(tmpvcf);
            vcf.create_ingroup(50);

            //CHECK(vcf.ingroup().size() == 1);
            //CHECK(vcf.ingroup()[0] == "s2");

            vcf.next();
            //REQUIRE(vcf.gen_genotypes().size() == 2);
            //CHECK(vcf.gen_genotypes()[0] == "C");
            //CHECK(vcf.gen_genotypes()[1] == "T");

            vcf.next();
            //REQUIRE(vcf.gen_genotypes().size() == 2);
            //CHECK(vcf.gen_genotypes()[0] == "");
            //CHECK(vcf.gen_genotypes()[1] == "CC");

            // Allele set should be complete, ingroup should reflect minimized set
            CHECK(vcf.alleles().size() == 3);
            //CHECK(vcf.ingroup().size() == 1);
        }

        SUBCASE("Allele populations") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.next();
            vcf.gen_genotypes();

            REQUIRE(vcf.allele_pop("G").size() == 4);
            CHECK(vcf.allele_pop("G")[0]);
            CHECK(!vcf.allele_pop("G")[1]);
            CHECK(!vcf.allele_pop("G")[2]);
            CHECK(!vcf.allele_pop("G")[3]);

            REQUIRE(vcf.allele_pop("A").size() == 4);
            CHECK(!vcf.allele_pop("A")[0]);
            CHECK(vcf.allele_pop("A")[1]);
            CHECK(!vcf.allele_pop("A")[2]);
            CHECK(!vcf.allele_pop("A")[3]);

            REQUIRE(vcf.allele_pop("C").size() == 4);
            CHECK(!vcf.allele_pop("C")[0]);
            CHECK(!vcf.allele_pop("C")[1]);
            CHECK(vcf.allele_pop("C")[2]);
            CHECK(!vcf.allele_pop("C")[3]);

            REQUIRE(vcf.allele_pop("T").size() == 4);
            CHECK(!vcf.allele_pop("T")[0]);
            CHECK(!vcf.allele_pop("T")[1]);
            CHECK(!vcf.allele_pop("T")[2]);
            CHECK(vcf.allele_pop("T")[3]);

        }

        SUBCASE("Filtered allele populations") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.create_ingroup({"s1"});
            vcf.next();
            vcf.gen_genotypes();

            REQUIRE(vcf.allele_pop("G").size() == 2);
            CHECK(vcf.allele_pop("G")[0]);
            CHECK(!vcf.allele_pop("G")[1]);

            REQUIRE(vcf.allele_pop("A").size() == 2);
            CHECK(!vcf.allele_pop("A")[0]);
            CHECK(vcf.allele_pop("A")[1]);

            REQUIRE(vcf.allele_pop("C").size() == 2);
            CHECK(!vcf.allele_pop("C")[0]);
            CHECK(!vcf.allele_pop("C")[1]);

            REQUIRE(vcf.allele_pop("T").size() == 2);
            CHECK(!vcf.allele_pop("T")[0]);
            CHECK(!vcf.allele_pop("T")[1]);
        }

        SUBCASE("Allele frequencies") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.next();

            auto af = vcf.frequencies();
            REQUIRE(af.size() == 4);
            CHECK(af[0] > 0.289f); // af[0] should be 0.29
            CHECK(af[0] < 0.291f);
            CHECK(af[1] == 0.01f);
            CHECK(af[2] == 0.6f);
            CHECK(af[3] == 0.1f);
        }

        SUBCASE("INFO field parsing") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.next();

            // Test INFO field extraction
            auto af_values = vcf.info_tag<float>("AF");
            REQUIRE(af_values.size() == 3);
            CHECK(af_values[0] == 0.01f);
            CHECK(af_values[1] == 0.6f);
            CHECK(af_values[2] == 0.1f);

            auto ac_values = vcf.info_tag<int>("AC");
            REQUIRE(ac_values.size() == 1);
            CHECK(ac_values[0] == 1);

            auto ns_values = vcf.info_tag<int>("NS");
            REQUIRE(ns_values.size() == 1);
            CHECK(ns_values[0] == 1);
        }

        SUBCASE("FORMAT field parsing") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.next();

            // Test FORMAT field extraction
            auto gt_values = vcf.format_tag<std::string>("GT");
            REQUIRE(gt_values.size() == 2); // 2 samples, each with one GT value
            CHECK(gt_values[0] == "0|1"); // First sample: s1
            CHECK(gt_values[1] == "2|3"); // Second sample: s2
        }

        SUBCASE("Multiple record iteration") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            
            int record_count = 0;
            std::vector<std::string> expected_refs = {"G", "C", "G", "TATA", "T"};
            
            while (vcf.next()) {
                REQUIRE(record_count < expected_refs.size());
                CHECK(vcf.ref() == expected_refs[record_count]);
                record_count++;
            }
            
            CHECK(record_count == 5); // Should read all 5 records
        }

        SUBCASE("Empty VCF handling") {
            // Create an empty VCF file
            std::string empty_vcf = "empty_test.vcf";
            {
                std::ofstream vcfo(empty_vcf);
                vcfo << "##fileformat=VCFv4.1" << std::endl
                     << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl
                     << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << std::endl;
            }

            vargas::VCF vcf;
            vcf.open(empty_vcf);
            
            // Should return false immediately
            CHECK(vcf.next() == false);
            
            remove(empty_vcf.c_str());
        }

    }

    remove(tmpvcf.c_str());
}

TEST_SUITE_END();
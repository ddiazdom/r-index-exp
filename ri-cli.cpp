//
// Created by Diaz, Diego on 12.11.2025.
//
#include "CLI11.hpp"
#include "internal/r_index.hpp"
#include <filesystem>
#include <chrono>
#include <random>

#define MEASURE(query, time_answer, query_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

#define MEASURE_VOID(query, time_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

std::vector<std::string> file2pat_list(std::string& pat_file, ulint &n_pats, ulint& pat_len){
    std::ifstream ifs(pat_file);
    std::string header;
    std::getline(ifs, header);
    n_pats = get_number_of_patterns(header);
    pat_len = get_patterns_length(header);
    std::cout<<"Searching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    std::vector<std::string> pat_list(n_pats);
    for(ulint i=0;i<n_pats;++i){
        pat_list[i].reserve(pat_len);
        for(ulint j=0;j<pat_len;++j){
            char c;
            ifs.get(c);
            pat_list[i].push_back(c);
        }
    }
    return pat_list;
}

// Helper: generate a random hex string for unique names
std::string random_hex(std::size_t length = 16) {
    static std::mt19937_64 rng{std::random_device{}()};
    static std::uniform_int_distribution<uint64_t> dist;

    std::ostringstream oss;
    while (oss.tellp() < static_cast<std::streamoff>(length)) {
        oss << std::hex << std::setw(16) << std::setfill('0') << dist(rng);
    }
    return oss.str().substr(0, length);
}

namespace fs = std::filesystem;
std::string create_tmp_dir(std::string prefix=""){
    fs::path temp_root;
    if (prefix.empty()) {
        temp_root = fs::temp_directory_path();
    }else {
        temp_root = prefix;
    }

    fs::path temp_dir = temp_root / fs::path("sri_"+random_hex(12));
    while (fs::exists(temp_dir)) {
        temp_dir = temp_root / fs::path("sri_"+random_hex(12));
    }
    fs::create_directories(temp_dir);

    return std::filesystem::absolute(temp_dir).lexically_normal().string();
}

struct arguments{
    std::string input_file;
    std::string output_file;
    std::string pat_file;
    std::string bigbwt_pref;
    size_t bytes_sa=5;
};

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

static void parse_app(CLI::App& app, struct arguments& args){

	auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(23);
    app.formatter(fmt);

    auto * build = app.add_subcommand("build");
    build->add_option("BIGBWT_PREFIX", args.bigbwt_pref, "BigBWT prefix containing the already computed BWT and SA samples")->required()->check(CLI::ExistingFile);
    build->add_option("-o,--output", args.output_file, "Output file where the index will be stored");

    auto * count = app.add_subcommand("count");
    count->add_option("INDEX", args.input_file, "Index file")->check(CLI::ExistingFile)->required();
    count->add_option("PAT_FILE", args.pat_file, "List of patterns")->check(CLI::ExistingFile)->required();

    auto * locate = app.add_subcommand("locate");
    locate->add_option("INDEX", args.input_file, "Index file")->check(CLI::ExistingFile)->required();
    locate->add_option("PAT_FILE", args.pat_file, "List of patterns")->check(CLI::ExistingFile)->required();

    app.require_subcommand(1,1);
}

template<class index_type>
void test_count(std::string input_file, std::string& pat_file){

    std::ifstream in(input_file);
    bool fast;
    in.read((char*)&fast,sizeof(fast));

    index_type index;
    index.load(in);

    uint64_t n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    std::pair<ulint, ulint> ans;
    size_t acc_count=0;
    for(auto const& p : pat_list) {
        MEASURE(index.count(p), acc_time, ans, std::chrono::nanoseconds)
        acc_count+=ans.second-ans.first;
    }

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"\tTotal number of occurrences "<<acc_count<<std::endl;
    std::cout<<"\t"<<double(acc_time)/double(n_pats)<<" nanosecs/pat"<<std::endl;
    std::cout<<"\t"<<double(acc_time)/double(acc_count)<<" nanosecs/occ"<<std::endl;
}

template<class index_type>
void test_locate(std::string input_file, std::string& pat_file){

    std::ifstream in(input_file);
    bool fast;
    in.read((char*)&fast,sizeof(fast));

    index_type index;
    index.load(in);

    uint64_t n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    size_t acc_count=0;
    for(auto const& p : pat_list) {
        std::vector<size_t> occ;
        MEASURE(index.locate_all(p), acc_time, occ, std::chrono::microseconds)
        acc_count+=occ.size();
    }

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"\tTotal number of occurrences "<<acc_count<<std::endl;
    std::cout<<"\t"<<double(acc_time)/double(n_pats)<<" microsecs/pat"<<std::endl;
    std::cout<<"\t"<<double(acc_time)/double(acc_count)<<" microsecs/occ"<<std::endl;
}

int main(int argc, char** argv) {

    arguments args;
    CLI::App app("Subsample r-index");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand("build")) {

        if(args.output_file.empty())  args.output_file = std::filesystem::path(args.bigbwt_pref).filename();
        if(fs::path(args.output_file).extension()!="ri") args.output_file += ".ri";

        std::cout<<"Building the r-index from the precomputed BWT/SA elements in "<<args.bigbwt_pref<<std::endl;
        auto idx = ri::r_index(args.bigbwt_pref);
        bool fast = false;//a parameter of the r-index, seems to be unused
        std::ofstream ofs(args.output_file);
        ofs.write((char*)&fast,sizeof(fast));

        idx.serialize(ofs);
        std::cout<<"The output index was stored in "<<args.output_file<<std::endl;
        ofs.close();
    } else if(app.got_subcommand("count")){
        test_count<ri::r_index<>>(args.input_file, args.pat_file);
    } else if(app.got_subcommand("locate")){
        test_locate<ri::r_index<>>(args.input_file, args.pat_file);
    } else {
        std::cerr<<" Unknown command "<<std::endl;
        exit(1);
    }
    return 0;
}
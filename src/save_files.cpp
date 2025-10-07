#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <cstdio>

// [[Rcpp::depends(RcppParallel, qs2)]]
// [[Rcpp::plugins(cpp14)]]

#include "qs2_external.h"

using namespace Rcpp;
using namespace RcppParallel;

namespace {

struct ParallelSaver : public Worker {
    const Rcpp::List& objects;
    const Rcpp::CharacterVector paths;
    const int compress_level;
    const bool shuffle;
    std::vector<std::string>& errors;
    const bool use_qdata;

    ParallelSaver(const Rcpp::List& objects_,
                  const Rcpp::CharacterVector& paths_,
                  int compress_level_, bool shuffle_,
                  std::vector<std::string>& errors_,
                  bool use_qdata_)
        : objects(objects_),
          paths(paths_),
          compress_level(compress_level_),
          shuffle(shuffle_),
          errors(errors_),
          use_qdata(use_qdata_) {
    }

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            try {
                SEXP obj = objects[i];

                size_t buffer_len = 0;
                unsigned char* buffer = nullptr;
                if (use_qdata) {
                    buffer = c_qd_serialize(obj, &buffer_len, compress_level, shuffle, false, 1);
                } else {
                    buffer = c_qs_serialize(obj, &buffer_len, compress_level, shuffle, 1);
                }

                if (buffer == nullptr) {
                    errors[i] = "Serialization returned NULL buffer";
                    continue;
                }

                Rcpp::String path_sexp = paths[i];
                const std::string path_str(path_sexp.get_cstring());

                FILE* handle = std::fopen(path_str.c_str(), "wb");
                if (handle == nullptr) {
                    if (use_qdata) {
                        c_qd_free(static_cast<void*>(buffer));
                    } else {
                        c_qs_free(static_cast<void*>(buffer));
                    }
                    errors[i] = "Failed to open file: " + path_str;
                    continue;
                }

                size_t written = std::fwrite(buffer, sizeof(unsigned char), buffer_len, handle);
                std::fclose(handle);
                if (use_qdata) {
                    c_qd_free(static_cast<void*>(buffer));
                } else {
                    c_qs_free(static_cast<void*>(buffer));
                }

                if (written != buffer_len) {
                    errors[i] = "Incomplete write: " + path_str;
                }
            } catch (std::exception& ex) {
                errors[i] = ex.what();
            } catch (...) {
                errors[i] = "Unknown error";
            }
        }
    }
};

} // namespace

//' Parallel qdata saver
//' 
//' @param objects List of data-only R objects to serialize (unsupported types become NULL, mirroring `qd_save`).
//' @param paths Character vector of output file paths (must match length of `objects`).
//' @param compress_level Compression level passed to qdata (default 3).
//' @param shuffle Whether to enable byte shuffling (default TRUE).
//' @param grain_size Chunk size for parallel processing (default 1).
//' @return NULL invisibly on success; character vector of error messages otherwise.
//' @export
// [[Rcpp::export]]
SEXP save_qd_cpp(const Rcpp::List& objects,
                 const Rcpp::CharacterVector& paths,
                 int compress_level = 3,
                 bool shuffle = true,
                 std::size_t grain_size = 1) {

    if (objects.size() != paths.size()) {
        Rcpp::stop("`objects` and `paths` must have the same length");
    }

    if (grain_size == 0) {
        Rcpp::stop("`grain_size` must be positive");
    }

    std::vector<std::string> errors(paths.size());
    ParallelSaver saver(
        objects,
        paths,
        compress_level,
        shuffle,
        errors,
        true
    );

    parallelFor(0, static_cast<std::size_t>(objects.size()), saver, grain_size);

    std::vector<std::string> reported;
    reported.reserve(errors.size());
    for (std::size_t i = 0; i < errors.size(); ++i) {
        if (!errors[i].empty()) {
            std::string msg = std::string(paths[i]) + " -> " + errors[i];
            Rcpp::Rcout << "save_qd_cpp: " << msg << '\n';
            reported.push_back(std::move(msg));
        }
    }

    if (!reported.empty()) {
        return Rcpp::wrap(reported);
    }

    return R_NilValue;
}



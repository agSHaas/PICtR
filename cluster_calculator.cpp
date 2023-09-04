#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector calculate_clusters(NumericMatrix cosine_similarities,
                                   List data_spleen_ly6c,
                                   int chunk_size) {
  
  int n_rows = cosine_similarities.nrow();
  CharacterVector result(n_rows);
  
  for (int row = 0; row < n_rows; ++row) {
    LogicalVector above_cutoff = cosine_similarities(row, _) > 0;
    
    std::vector<std::string> cells;
    for (int i = 0; i < above_cutoff.size(); ++i) {
      if (above_cutoff[i]) {
        cells.push_back(data_spleen_ly6c[i]);
      }
    }
    
    std::unordered_map<std::string, int> cluster_count;
    std::string max_cluster;
    int max_count = 0;
    
    for (const auto& cell : cells) {
      std::string cluster = data_spleen_ly6c[cell];
      cluster_count[cluster]++;
      if (cluster_count[cluster] > max_count) {
        max_cluster = cluster;
        max_count = cluster_count[cluster];
      }
    }
    
    result[row] = max_cluster;
  }
  
  return result;
}


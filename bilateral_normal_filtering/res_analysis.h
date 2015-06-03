#ifndef RES_ANALYSIS_H
#define RES_ANALYSIS_H

#include <map>
#include <zjucad/matrix/matrix.h>

namespace zsw
{
  void analysisFaceNormalEnergy(const std::string &vtk_file, const std::string &output_vtk);
  void calFaceNormalEnergy(const std::multimap<size_t, size_t> &one_ring_f, const zjucad::matrix::matrix<double> &fnormal,
                            zjucad::matrix::matrix<double> &energy);
}

#endif /* RES_ANALYSIS_H */

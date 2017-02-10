#include <algorithm>
#ifdef __cplusplus
extern "C"
{
#endif

  void nthelement(double *x, int num, int n)
  {
    std::nth_element(x, x+num, x+n);
  }

#ifdef __cplusplus
}
#endif


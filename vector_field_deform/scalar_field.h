#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

namespace zsw {

  class Function
  {
  public:
    virtual double val(const double *x) = 0;
    virtual void jac(const double *x, double *g) = 0;
  };

  /**
   * function e(x) and f(x)
   */
  class LinearScalarField  final : public Function
  {
  public:
    LinearScalarField(const double *u, const double *c);
    virtual double val(const double *x);
    virtual void jac(const double *x, double *g);
    ~LinearScalarField();
  private:
    double u_[3];
    double c_[3];
  };

  /**
   * f(x) = [a.cross(x-c)]^2
   */
  class QuadraticScalarField final : public Function
  {
  public:
    QuadraticScalarField(const double *a, const double *c);
    double val(const double *x);
    void jac(const double *x, double *g);
  private:
    double a_[3];
    double c_[3];
  };

  /**
   * e(x) = [a.dot(x-c)]^2
   */
  class QuadraticScalarField2 final : public Function
  {
  public:
    QuadraticScalarField2(const double *a, const double *center);
    double val(const double *x);
    void jac(const double *x, double *g);
  private:
    double a_[3];
    double center_[3];
  };

  /**
   * function b(r)
   */
  class BlendFunc final : public Function
  {
  public:
    BlendFunc(const double ri, const double ro);
    double val(const double *r);
    void jac(const double *r, double *g);
    ~BlendFunc();
  private:
    double ri_;
    double ro_;
  };

  /**
   * function r(x)
   */
  class RegionFunc
  {
  public:
    enum REGION_TYPE {
      INNER_REGION,
      BLENDER_REGION,
      OUTER_REGION
    };
    virtual double val(const double*x) = 0;
    virtual void jac(const double *x, double *g) = 0;
    virtual REGION_TYPE judgeRegion(const double *x) = 0;
    double getRi() const { return ri_; }
    double getRo() const { return ro_; }
  protected:
    double ri_;
    double ro_;
  };

  class SphereRegionFunc final : public RegionFunc
  {
  public:
    SphereRegionFunc(const double ri, const double ro, const double *center);
    double val(const double*x);
    void jac(const double *x, double *g);
    REGION_TYPE judgeRegion(const double *x);
    const double * getCenter() const { return center_; }
  private:
    double center_[3];
  };

  class IsosurfacesRegionFunc final : public RegionFunc
  {
  public:
    IsosurfacesRegionFunc(const double *b, const double *center, const double ri, const double ro);
    double val(const double*x);
    void jac(const double *x, double *g);
    REGION_TYPE judgeRegion(const double *x);
  private:
    double center_[3];
    double b_[3];
  };
}

#endif /* SCALAR_FIELD_H */

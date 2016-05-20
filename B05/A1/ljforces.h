#ifndef ljforces_H
#define ljforces_H
class LJForces {
    public:
    Eigen::MatrixXd kraft(Eigen::MatrixXd particleinfo, double L);
};
#endif

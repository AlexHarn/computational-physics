#ifndef periodRB_H
#define periodRB_H
class PeriodRB {
    public:
    Eigen::Vector2d umklapp(Eigen::Vector2d tempRB, double L);
    Eigen::Vector2d kurzerWeg(Eigen::Vector2d dr, double L);
};
#endif

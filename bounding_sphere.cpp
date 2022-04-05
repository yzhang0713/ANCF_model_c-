#include "bounding_sphere.h"

int check_disjoint(const bounding_sphere & bs1, const bounding_sphere & bs2) {
    // Get distance between bounding sphere centers
    double distance = (bs1.center_point - bs2.center_point).norm();
    return distance > bs1.radius + bs2.radius;
}

int bounding_sphere::is_bs_disjoint(const bounding_sphere & other) {
    return check_disjoint(*this, other);
}
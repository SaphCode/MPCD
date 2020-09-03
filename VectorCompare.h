#pragma once

#ifndef VECTORCOMPARE_H
#define VECTORCOMPARE_H

#include <functional>
#include <Eigen/Dense>
#include <map>

namespace std
{
    template <> struct less<Eigen::Vector2i>
    {
        bool operator() (const Eigen::Vector2i& lhs, const Eigen::Vector2i& rhs) const
        {
            assert(lhs.size() == rhs.size());
            for (int i = 0; i < lhs.size(); i++) {
                if (lhs[i] < rhs[i]) return true;
                if (lhs[i] > rhs[i]) return false;
            }
            return false;
        }
    };
}

namespace MPCD {
    typedef std::map<Eigen::Vector2i, Eigen::Vector2d, std::less<Eigen::Vector2i>, Eigen::aligned_allocator < std::pair<const Eigen::Vector2i, Eigen::Vector2d>>> vectorMap;
}
#endif
#pragma once
// #include <omp.h>
#include "Eigen/Dense"

class Photon {
public:
    Eigen::Vector3f pos;                    // photon position
    Eigen::Vector3f dir;                    // incident direction
    int axis;                               // splitting axis
public:
    Photon();                               // constructor

    friend class Map;
};

class Nearest_photons {
public:
    int max_num;                            // number of nearest photons required
    int curr_num;                           // number of nearest photons found
    bool built;                             // whether the heap has been built up
    Eigen::Vector3f pos;                    // position to search photons around
    float* dist;                            // array of distance from a specific photon to the required position
    Photon** photons;                       // array of photons found
public:
    Nearest_photons(                        // constructor
        int n,
        Eigen::Vector3f p,
        float d);
    ~Nearest_photons();                     // destructor
    Photon** get_photons() const;           // photons retriever

    friend class Map;
};

class Map {
public:
    Photon* photons;                        // array of photons
    int stored_photons;                     // number of current photons
    int max_photons;                        // photons capacity
    Eigen::Vector3f light_power;            // radiance of the light (to calculate photon power)
    Eigen::Vector3f bbox_min;               // smallest coordinate of all photon position
    Eigen::Vector3f bbox_max;               // largest coordinate of all photon position
    void balance_segment(                   // balance the array (current root at root) from start to end
        Photon** out,
        Photon** in,
        int root,
        int start,
        int end);
public:
    Map(                                                                    // constructor
        int max_photons,
        Eigen::Vector3f light_power);
    ~Map();                                                                 // destructor
    void store(                                                             // call to store photons to photons array
        const Eigen::Vector3f& pos,
        const Eigen::Vector3f& dir);
    static Eigen::Vector3f photon_dir(const Photon* p);                   // direction retriever
    void balance();                                                         // call to build kd-tree from a flat array
    void locate_photons(                                                    // k-nearest neighbor algorithm
        Nearest_photons* np,
        int root = 1);
};

Photon::Photon() :
    pos(Eigen::Vector3f::Zero()),
    dir(Eigen::Vector3f::Zero()),
    axis(-1) {}

Nearest_photons::Nearest_photons(
    int n,
    Eigen::Vector3f p,
    float d) :
    max_num(n),
    curr_num(0),
    built(false),
    pos(std::move(p)) {
    dist = new float[n + 1];                  // allocate memories for the heap
    photons = new Photon * [n + 1];             // allocate memories for the heap
    dist[0] = d * d;
}

Nearest_photons::~Nearest_photons() {
    delete[] dist;
    delete[] photons;
}

Photon** Nearest_photons::get_photons() const {
    return photons;
}

Map::Map(
    int max_photons,
    Eigen::Vector3f light_power) :
    stored_photons(0),
    max_photons(max_photons),
    light_power(std::move(light_power)) {
    photons = new Photon[max_photons + 1];
    for (int i = 0; i < 3; i++) {
        bbox_min[i] = std::numeric_limits<float>::max();
        bbox_max[i] = -1 * std::numeric_limits<float>::max();
    }
}

Map::~Map() {
    delete[] photons;
}

void Map::store(
    const Eigen::Vector3f& pos,
    const Eigen::Vector3f& dir) {
    if (stored_photons == max_photons)                  // array is already full
        return;

    stored_photons++;                                  // add a new photon
    Photon* p = &photons[stored_photons];               // retrieve the back position

    p->pos = pos;                                       // set photon position
    p->dir = dir.normalized();                          // set incident direction

    bbox_min = bbox_min.cwiseMin(pos);                  // enlarge the bounding box lower bound
    bbox_max = bbox_max.cwiseMax(pos);                  // enlarge the bounding box upper bound
}

Eigen::Vector3f Map::photon_dir(const Photon* p) {
    return p->dir;
}

void Map::balance_segment(
    Photon** out,
    Photon** in,
    int root,
    int start,
    int end) {
    int axis;
    Eigen::Vector3f diff = bbox_max - bbox_min;                                 // calculate bounding box difference
    float max = std::max(std::max(diff.x(), diff.y()), diff.z());                 // find out axis of max difference
    if (diff.x() == max)
        axis = 0;
    else if (diff.y() == max)
        axis = 1;
    else
        axis = 2;

    int left = start;                                                           // calculate k-value for compact left balancing binary search trees
    int right = end;
    int median = 1;
    while (4 * median <= end - start + 1)
        median *= 2;
    if (3 * median <= end - start + 1) {
        median *= 2;
        median += start - 1;
    }
    else
        median = end - median + 1;

    while (right > left) {                                                      // quick select algorithm
        float v = in[right]->pos[axis];
        int i = left;                                                           // index of left guard
        int j = right - 1;                                                      // index of right guard
        while (true) {
            while (in[i]->pos[axis] <= v && i < right)                          // ensure all element to the left of median is smaller
                i++;                                                           // left guard move to right
            while (in[j]->pos[axis] >= v && j > left)                           // ensure all element to the right of median is larger
                j--;                                                           // right guard move to left
            if (i >= j)                                                         // break when sorted
                break;
            std::swap(in[i], in[j]);                                       // swap the first two unsorted element
        }
        std::swap(in[i], in[right]);                                       // put standard element to the right place
        if (i > median)                                                         // reduce search field to the right half
            right = i - 1;
        else                                                                    // reduce search field to the left half
            left = i + 1;
    }

    out[root] = in[median];                                                     // set elements for current root
    out[root]->axis = axis;

    if (median > start) {                                                       // if any photons in left sub tree
        if (start < median - 1) {
            float tmp = bbox_max[axis];
            bbox_max[axis] = out[root]->pos[axis];                              // reduce bounding box
            balance_segment(out, in, 2 * root, start, median - 1);             // balance left sub tree
            bbox_max[axis] = tmp;                                               // recover bounding box for the next call
        }
        else {                                                                // if only one photon in left sub tree
            out[2 * root] = in[start];
        }
    }

    if (median < end) {                                                         // if any photons in right sub tree
        if (median + 1 < end) {
            float tmp = bbox_min[axis];
            bbox_min[axis] = out[root]->pos[axis];                              // reduce bounding box
            balance_segment(out, in, 2 * root + 1, median + 1, end);            // balance right sub tree
            bbox_min[axis] = tmp;                                               // recover bounding box
        }
        else {                                                                // if only on photon in right sub tree
            out[2 * root + 1] = in[end];
        }
    }
}

void Map::balance() {
    if (stored_photons > 1) {
        auto** tmp1 = new Photon * [stored_photons + 1];
        auto** tmp2 = new Photon * [stored_photons + 1];

#pragma omp parallel for
        for (int i = 0; i <= stored_photons; i++)
            tmp2[i] = &photons[i];

        balance_segment(tmp1, tmp2, 1, 1, stored_photons);

        delete[] tmp2;
        auto* tmp3 = new Photon[stored_photons + 1];

#pragma omp parallel for
        for (int i = 1; i <= stored_photons; i++)
            tmp3[i] = *tmp1[i];

        delete[] tmp1;
        delete[] photons;
        photons = tmp3;
    }
}

void Map::locate_photons(
    Nearest_photons* np,
    int root) {
    Photon* p = &photons[root];
    if (root < stored_photons / 2 - 1) {                                                                        // if current node is not leaf node
        float dist_to_bound = np->pos[p->axis] - p->pos[p->axis];                                           // calculate vertical distance to boundary
        if (dist_to_bound > 0.0f) {                                                                         // if position required is in the right half
            locate_photons(np, 2 * root + 1);                                                               // call for the right sub tree
            if (dist_to_bound * dist_to_bound < np->dist[0]) {                                                // call for the left sub tree if necessary
                locate_photons(np, 2 * root);
            }
        }
        else {                                                                                            // if position required is in the left half
            locate_photons(np, 2 * root);                                                                 // call for the left sub tree
            if (dist_to_bound * dist_to_bound < np->dist[0]) {                                                // call for the right sub tree if necessary
                locate_photons(np, 2 * root + 1);
            }
        }
    }

    float dist_to_photon = (p->pos - np->pos).norm();                                                       // calculate distance to the photon
    if (dist_to_photon < np->dist[0]) {                                                                     // if condition satisfied
        if (np->curr_num < np->max_num) {                                                                   // when heap is not full
            np->curr_num++;                                                                                // add a new photon
            np->dist[np->curr_num] = dist_to_photon;                                                        // add element distance
            np->photons[np->curr_num] = p;                                                                  // add element photon
        }
        else {
            int parent, child;
            if (!np->built) {                                                                               // if heap is not set up
                float dist;
                Photon* ptr;
                for (int i = np->max_num / 2; i > 0; i--) {                                                  // Floyd algorithm
                    parent = i;
                    dist = np->dist[parent];
                    ptr = np->photons[parent];
                    while (parent <= np->max_num / 2) {
                        child = 2 * parent;
                        if (np->dist[child] < np->dist[child + 1] && child + 1 <= np->max_num)                  // if right child is larger
                            child++;
                        if (dist >= np->dist[child])                                                        // parent larger than both children, no need to modify
                            break;
                        np->dist[parent] = np->dist[child];                                                 // exchange parent and the larger child if necessary
                        np->photons[parent] = np->photons[child];
                        parent = child;
                    }
                    np->dist[parent] = dist;
                    np->photons[parent] = ptr;
                }
                np->built = true;                                                                           // heap is built, set the flag to true
            }
            parent = 1;
            for (child = 2; child <= np->max_num; child *= 2) {
                if (np->dist[child] < np->dist[child + 1] && child != np->max_num)                            // if right child is larger
                    child++;
                if (dist_to_photon > np->dist[child])                                                       // larger than both children, no need to step through
                    break;
                np->dist[parent] = np->dist[child];                                                         // exchange parent node and the larger child node
                np->photons[parent] = np->photons[child];
                parent = child;
            }
            if (dist_to_photon < np->dist[parent]) {
                np->photons[parent] = p;
                np->dist[parent] = dist_to_photon;
            }
            np->dist[0] = np->dist[1];
        }
    }
}
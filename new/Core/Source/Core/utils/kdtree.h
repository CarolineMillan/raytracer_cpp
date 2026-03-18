#pragma once

#include <vector>
#include <queue>
#include <algorithm>
#include <limits>
#include <cmath>

#include "photon.h"   // defines class Photon { Vertex position; ... }
#include "vector.h"       // defines Vector

// KD-tree node storing a Photon pointer
struct KDNode {
    Photon* photon;
    int     axis;   // splitting axis: 0=x, 1=y, 2=z
    KDNode* left;
    KDNode* right;
};

class KDTree {
public:
    // Build a KD-tree from a list of Photon* pointers
    explicit KDTree(std::vector<Photon*>& photons) {
        root = buildRecursive(photons, 0, static_cast<int>(photons.size()), 0);
    }

    // Find the k nearest photons to 'query' and store in 'out'
    void kNearest(const Vertex& query, int k, std::vector<Photon*>& out) const {
        // max-heap of (squared distance, Photon*)
        std::priority_queue<std::pair<float, Photon*>> best;
        knnRecursive(root, query, k, best);
        out.clear();
        while (!best.empty()) {
            out.push_back(best.top().second);
            best.pop();
        }
        std::reverse(out.begin(), out.end());  // nearest-first
    }

    ~KDTree() {
        freeRecursive(root);
    }

private:
    KDNode* root = nullptr;

    // Recursive build: [begin, end) range in photons, depth for axis
    KDNode* buildRecursive(std::vector<Photon*>& pts, int begin, int end, int depth) {
        if (begin >= end) return nullptr;
        int axis = depth % 3;
        int mid = (begin + end) / 2;
        auto cmp = [axis](Photon* a, Photon* b) {
            switch (axis) {
                case 0: return a->position.x < b->position.x;
                case 1: return a->position.y < b->position.y;
                default: return a->position.z < b->position.z;
            }
        };
        std::nth_element(pts.begin() + begin, pts.begin() + mid, pts.begin() + end, cmp);
        KDNode* node = new KDNode{ pts[mid], axis, nullptr, nullptr };
        node->left  = buildRecursive(pts, begin, mid, depth + 1);
        node->right = buildRecursive(pts, mid + 1, end, depth + 1);
        return node;
    }

    // Recursive k-NN search
    void knnRecursive(KDNode* node,
                      const Vertex& query,
                      int k,
                      std::priority_queue<std::pair<float, Photon*>>& best) const {
        if (!node) return;
        // compute squared distance to this photon
        float dx = query.x - node->photon->position.x;
        float dy = query.y - node->photon->position.y;
        float dz = query.z - node->photon->position.z;
        float dist2 = dx*dx + dy*dy + dz*dz;
        if ((int)best.size() < k) {
            best.emplace(dist2, node->photon);
        } else if (dist2 < best.top().first) {
            best.pop();
            best.emplace(dist2, node->photon);
        }
        // decide which subtree to search first
        float diff;
        switch (node->axis) {
            case 0: diff = query.x - node->photon->position.x; break;
            case 1: diff = query.y - node->photon->position.y; break;
            default: diff = query.z - node->photon->position.z; break;
        }
        KDNode* first  = (diff < 0) ? node->left : node->right;
        KDNode* second = (diff < 0) ? node->right : node->left;
        knnRecursive(first, query, k, best);
        // if hypersphere crosses splitting plane, check other side
        if ((int)best.size() < k || diff*diff < best.top().first) {
            knnRecursive(second, query, k, best);
        }
    }

    // Recursively free all nodes
    void freeRecursive(KDNode* node) {
        if (!node) return;
        freeRecursive(node->left);
        freeRecursive(node->right);
        delete node;
    }
};

#pragma once

#include <parallel/algorithm>
#include "Node.h"
#include "wtime.h"


#if 0
#define SLOW
#endif

struct Tree
{
  typedef boundary<float> Boundary;

  Particle::Vector tree;
  Boundary BBox;  /* bounding box */
  std::vector<Node *> leafArray;

  struct cmp_particle_key { bool operator() (const Particle &a, const Particle &b) {return a.key.val < b.key.val;} };


  Tree(const Particle::Vector &ptcl_in)
  {
    const double t0 = wtime();

    std::vector<Particle> &ptcl = Node::ptcl;
    ptcl = ptcl_in;
    const int nbody = ptcl_in.size();

    /* import particles and compute the Bounding Box */
    for (int i = 0; i < nbody; i++)
      BBox.merge(Boundary(ptcl[i].pos));

    std::cerr << BBox.min << std::endl;
    std::cerr << BBox.max << std::endl;

    const vec3  vsize = BBox.hlen();
    const float rsize = std::max(vsize.x, std::max(vsize.x, vsize.y)) * 2.0f;

    float rsize2 = 1.0;
    while (rsize2 > rsize) rsize2 *= 0.5;
    while (rsize2 < rsize) rsize2 *= 2.0;

    /* now build the tree */

    for (int i = 0; i < nbody; i++)
      ptcl[i].compute_key(BBox.min, rsize2);

    __gnu_parallel::sort(ptcl.begin(), ptcl.end(), cmp_particle_key());

    Node::Node_heap.push_back(Node());
    Node &root = Node::Node_heap[0];
    root.size = rsize2;
    for (int i = 0; i < nbody; i++)
      root.push_particle(i, 60);

    root.make_boundary();

    root.find_group_Node(Node::NLEAF, leafArray);

    const double t1 = wtime();
    fprintf(stderr, " -- Tree build is done in %g sec [ %g ptcl/sec ]\n",  t1 - t0, nbody/(t1 - t0));

  };
};



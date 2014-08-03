#include <iostream>
#include <omp.h>
#include "Particle.h"
#include "boundary.h"

struct Node
{
  static const int NLEAF = 64;
  static std::vector<Particle> ptcl;
  static std::vector<Node> Node_heap;
  static std::vector<std::pair<Node*, Node*> > pair_list;
  static void allocate(int nptcl, int nNode){
    ptcl.reserve(nptcl);
    Node_heap.reserve(nNode);
  }
  static void clear(){
    ptcl.clear();
    Node_heap.clear();
  }
  typedef boundary<float> Boundary;

  int np;     // number of Particle;
  int depth;
  float size;
  int pfirst; // first Particle
  int cfirst; // first child
  Boundary bound_inner;
  Boundary bound_outer;

  Node()           : np(0), depth(     0), pfirst(-1), cfirst(-1) {}
  Node(int _depth, float _size) : np(0), depth(_depth), size(_size), pfirst(-1), cfirst(-1) {}

  bool is_leaf() const
  {
    return np < NLEAF;
  }
  void push_particle(const int paddr, const int rshift)
  {
    assert(rshift >= 0);
    if(!is_leaf())
    { // assign recursively
      int ic = ptcl[paddr].octkey(rshift);
      Node &child = Node_heap[cfirst + ic];
      child.push_particle(paddr, rshift-3);
    }
    else
    {
      if(-1 == pfirst)
      {
        assert(0 == np);
        pfirst = paddr;
      }
    }
    np++;
    if(np == NLEAF)
    { // shi's just become a mother
      assert(pfirst >= 0);
      cfirst = Node_heap.size();
#if 0
      for(int ic=0; ic<8; ic++){
        Node_heap.push_back(Node(1+depth));
      }
#else
      size_t new_size = Node_heap.size() + 8;
      assert(Node_heap.capacity() >= new_size);
      Node_heap.resize(new_size, Node(1+depth,size/2.0));
#endif
      for(int addr = pfirst; addr < pfirst+np; addr++)
      {
        int ic = ptcl[addr].octkey(rshift);
        Node &child = Node_heap[cfirst + ic];
        child.push_particle(addr, rshift-3);
      }
    }
  }
  void dump_tree(
      int level,
      std::ostream &ofs = std::cout) const{
    if(is_leaf()){
      for(int ip=0; ip<np; ip++){
        const Particle &p = ptcl[ip+pfirst];
        for(int i=0; i<level; i++) ofs << " ";
        ofs << p.pos << std::endl;
      }
      ofs << std::endl;
    }else{
      for(int i=0; i<level; i++) ofs << ">";
      ofs << std::endl;
      for(int ic=0; ic<8; ic++){
        const Node &child = Node_heap[cfirst + ic];
        child.dump_tree(level+1, ofs);
      }
    }
  }
  void make_boundary(){
    if(is_leaf()){
      for(int ip=0; ip<np; ip++){
        const Particle &p = ptcl[ip+pfirst];
        bound_inner.merge(Boundary(p.pos));
        bound_outer.merge(Boundary(p.pos, p.get_h()));
      }
    }else{
      for(int ic=0; ic<8; ic++){
        Node &child = Node_heap[cfirst + ic];
        if(child.np > 0){
          child.make_boundary();
          bound_inner.merge(child.bound_inner);
          bound_outer.merge(child.bound_outer);
        }
      }
    }
  }
  
  void find_group_Node(
      int ncrit,
      std::vector<Node *> &group_list){
    if (np == 0)
      return;
    if(np < ncrit){
      group_list.push_back(this);
    }else{
      for(int ic=0; ic<8; ic++){
        Node &child = Node_heap[cfirst + ic];
        child.find_group_Node(ncrit, group_list);
      }
    }
  }

};

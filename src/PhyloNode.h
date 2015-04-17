
#ifndef GUARD_PhyloNode
#define GUARD_PhyloNode


#include "global.hpp"

#include "tools.hpp"


typedef int32_t tax_id;
typedef int8_t rank_t;
#define MAX_TREE_DEPTH 256





class PhyloNode;
class Match_hit;


class PhyloNode {
	public:

	tax_id ncbi_tax_id;

	std::vector<PhyloNode * > * children;
	std::vector<Match_hit * > * hits;




	int positives;
	int non_positives;


	double f_measure;




	PhyloNode(tax_id ncbi_tax_id)
		: ncbi_tax_id(ncbi_tax_id), children(0), hits(0), positives(0), non_positives(0) {};

	~PhyloNode();

	PhyloNode * getNode(tax_id ncbi_tax_id); // adds node if not exists, returns node.
	PhyloNode * findNode(tax_id ncbi_tax_id); // returns node if exists


	void addHit(Match_hit * hit);
	tax_id * getLCA();
	tax_id * getLCA(tax_id * lca, int pos);

	void printTree(std::ofstream& stream);
	void setHitsZero();
	void computeFMeasure(int positives_total, int best_length, std::list<PhyloNode * > * good_f_measure_nodes);

};


#endif

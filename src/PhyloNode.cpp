

#include "PhyloNode.h"



using namespace std;

void PhyloNode::printTree(ofstream& stream) {
	stream << "(";



	if (this->children != 0) {
		vector<PhyloNode * >::iterator it;

		bool first = true;
		for (it = this->children->begin(); it < this->children->end(); it++) {
			if (first == false) stream << ",";

			(*it)->printTree(stream);

			first = false;
		}
	}

	if (this->hits != 0) {
		if (this->children != 0) stream << ",";

		vector<Match_hit * >::iterator it;

		bool first = true;
		for (it = this->hits->begin(); it < this->hits->end(); it++) {
			if (first == false) stream << ",";

//			if ((*it)->dist == -1) {
//				stream << "real species";
//			} else {
//
//				if ((*it)->strongness >= 10) {
//					int d = (*it)->dist*100;
//					for (int i = 0; i < d ; i++) stream << 'x';
//				} else {
//					stream << "???";
//				}
//			}

			first = false;
		}

	}

	stream << ")"<< this->ncbi_tax_id;
}

void PhyloNode::setHitsZero() {

	if (this->children != 0) {
		vector<PhyloNode * >::iterator it;

		for (it = this->children->begin(); it < this->children->end(); it++) {
			(*it)->setHitsZero();
		}
	}


	if (this->hits != 0) {
		delete this->hits;
		this->hits = 0;
	}

}




PhyloNode::~PhyloNode(){
	if (this->children != 0) {
		DeleteContainerWithPointers(this->children);
	}

	if (this->hits != 0) {
		DeleteContainerWithPointers(this->hits);
	}
}


tax_id * PhyloNode::getLCA() {
	tax_id * lca = new tax_id[MAX_TREE_DEPTH];

	return this->getLCA(lca, 0);
}

tax_id * PhyloNode::getLCA(tax_id * lca, int pos) {

	lca[pos] = this->ncbi_tax_id;
	lca[pos+1] = -1;

	if (this->hits != 0) {
			return lca;
	}

	if (this->children == 0) {
			cerr << "no children !?!" << endl;
			exit(1);
	}

	if (this->children->size() > 1) {
		// here the tree branches !
		return lca;
	} else {
		// no branching, go one node deeper.
		return this->children->front()->getLCA(lca, pos+1);
	}
}

PhyloNode * PhyloNode::getNode(tax_id ncbi_tax_id) {

	PhyloNode * pn = 0;

	if (this->children == 0) {
		this->children = new vector<PhyloNode * >;
	} else {
		pn = findNode(ncbi_tax_id);

		if (pn != 0) {
			return pn;
		}

	}

	pn = new PhyloNode(ncbi_tax_id);
	this->children->push_back(pn);

	return pn;
}

PhyloNode * PhyloNode::findNode(tax_id ncbi_tax_id) {

	if (this->children == 0) return 0;
	vector<PhyloNode * >::iterator it;

	for (it = this->children->begin(); it < this->children->end(); it++) {
			if ((*it)->ncbi_tax_id == ncbi_tax_id) {
					return *it;
			}
	}

	return 0;
}

void PhyloNode::addHit(Match_hit * hit){
	if (this->hits == 0) {
			this->hits = new vector<Match_hit * >;
	}

	this->hits->push_back(hit);
}

void PhyloNode::computeFMeasure(int positives_total, int best_length, std::list<PhyloNode * > * good_f_measure_nodes) {



	if (this->children != 0) {
		vector<PhyloNode * >::iterator it;

		for (it = this->children->begin(); it < this->children->end(); it++) {

			(*it)->computeFMeasure(positives_total, best_length, good_f_measure_nodes);
			//this->positives+=(*it)->positives;
			//this->non_positives+=(*it)->non_positives;
		}
	}

	if (this->hits != 0) {
		vector<Match_hit * >::iterator it;

		for (it = this->hits->begin(); it < this->hits->end(); it++) {
			//this->positives+=(*it)->positives;
			//this->non_positives+= (best_length - (*it)->positives);
		}

	}

	// precision= TP/P
	// recall= TP/(TP+FN)

	const double F_MEASURE_BETA=2; // = 2 ---> F2 measure weights recall twice as much as precision

	double beta_sqr = F_MEASURE_BETA*F_MEASURE_BETA;


//double precision = this->positives / (double) (this->positives + this->non_positives);
	double precision;
//double recall =  this->positives / (double) positives_total;
double recall;


	this->f_measure= (1+beta_sqr)* ((precision*recall)/((beta_sqr*precision)+recall));

//	cout << "------------------- "<< endl;
//	cout << "ncbi_tax_id: " << this->ncbi_tax_id << endl;
//	cout << "precision: " << precision << endl;
//	cout << "recall: " << recall << endl;
//	cout << "f_measure: " << this->f_measure << endl;
//



//cout << "XXXXXXXXXXXXXXXXXXXXXXXXx: " << good_f_measure_nodes->size() << endl;
	if (good_f_measure_nodes->size() < 10) {

		list<PhyloNode * >::iterator it;
		it = good_f_measure_nodes->begin();
		while (1) {
			if (it == good_f_measure_nodes->end()) {
				good_f_measure_nodes->push_back(this);
				break;
			}

			if (this->f_measure > (*it)->f_measure) {
				good_f_measure_nodes->insert(it, this);
				break;
			}

			it++;
		}



	} else {

		if (this->f_measure > good_f_measure_nodes->back()->f_measure) {
			list<PhyloNode * >::iterator it;
			it = good_f_measure_nodes->begin();

			while (this->f_measure <= (*it)->f_measure) {
					it++;
			}
			good_f_measure_nodes->insert(it, this);

			// delete last element (with lowest f-measure):
			good_f_measure_nodes->pop_back();
		}

	}



}

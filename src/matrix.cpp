

#include "matrix.hpp"


Class_Substitution_Matrix::~Class_Substitution_Matrix(){
	int size = aminoacid_alphabetSize+1;
	for (int i=0; i<size; i++) {
			delete [] this->matrix[i];
	}
	delete [] this->matrix;

	//delete [] this->self_scores;
	//delete [] this->counts;
	//delete [] this->neighbors;
	//delete [] this->neighbors_scores;
}





//Class_Substitution_Matrix::Class_Substitution_Matrix() {
//
//}

Class_Substitution_Matrix::Class_Substitution_Matrix(string * filename){

	aminoacid_alphabetSize=20; // 0-19
	gapChar = aminoacid_alphabetSize;
	aminoacid_ascii2num = new int [256];

	aminoacid_num2ascii = new char[aminoacid_alphabetSize+1];
	for (int i=0; i<= aminoacid_alphabetSize; i++) {aminoacid_num2ascii[i] = '-';}  // all upper case

	// hmm: ACDEFGHIKLMNPQRSTVWY
	for (int i=0; i< 256; i++) {aminoacid_ascii2num[i] = gapChar;}	  // initialize with invalid/gap character
	// no amino acids: O,U,B,J,X,Z - treat as gap...
	int j=0;
	for (int i='A'; i<= 'Z'; i++) {
		if (i=='B' || i=='J' || i=='O' || i=='U' || i=='X' || i=='Z') {
			aminoacid_ascii2num[i] = gapChar;
		} else {
			aminoacid_ascii2num[i] = j;
			aminoacid_num2ascii[j] = i;

			j++;
		}
	}



	this->matrix = readMatrix(filename);

}

int Class_Substitution_Matrix::getScore(char a, char b) {

	if (a == '?' || b == '?') {
		return 0;
	}

	int a_int = aminoacid_ascii2num[a];
	int b_int = aminoacid_ascii2num[b];

	if (a_int == gapChar && b_int == gapChar) {
		return 0;
	}

	return matrix[a_int][b_int];
}

int ** Class_Substitution_Matrix::readMatrix(string * filename){
	string line;
	string * alphabet=0;
	int size = aminoacid_alphabetSize+1;
	int ** substitution_matrix = new int*[size];

	for (int i=0; i<size; i++) {
			substitution_matrix[i]=new int[size];
			for (int j=0; j<size; j++) {
				substitution_matrix[i][j]=INT_MAX;
			}
	}

	int score;
	ifstream matrixStream(filename->c_str());

	if (! matrixStream.is_open())
	{
		cerr << "Error reading file " << *filename << endl;
		exit(1);
	}

	while (getline(matrixStream, line)) {

		if (line[0] == '#') {
			continue;
		}

		if (alphabet == 0) {
				alphabet=new string(line);
				int i=0;
				while ((i = alphabet->find(' ', i)) != -1) {
					alphabet->erase(i,1);
				}
			continue;
		}


		if ((ulong) aminoacid_ascii2num[(uchar) line[0]] == gapChar && (line[0]!='*')) {
			// we don't support index amino acid.
			continue;
		}

		// process matrix line

		int i=1;
		int j=1;
		int column=0;

		while (j < (int) line.length()-1) {

			//search next number
			while (line[i] == ' ') {
				i++;
			}

			j=line.find(' ',i);
			if (j==-1) {	// last entry in line
				j=line.length()-1;
			}




			if ((ulong) aminoacid_ascii2num[(uchar) alphabet->at(column)] == gapChar && (alphabet->at(column) != '*')) {
				// we don't support this amino acid.
				column++;
				i=j;
				continue;
			}
			string sub = line.substr(i-1,j-i+1);
			score = str2int(sub);
				//cout << column << " (" << score << ")" << endl;
			int x = aminoacid_ascii2num[(uchar) line[0]];
			int y = aminoacid_ascii2num[(uchar) alphabet->at(column)];
			substitution_matrix[x][y]=score;

			column++;
			i=j;
		} // end while

	}// end while


		//cout << *alphabet << endl;
		//exit(0);
		//str2int()


	matrixStream.close();
	delete alphabet;

	checkMatrix(substitution_matrix);
	return substitution_matrix;
}

void Class_Substitution_Matrix::printMatrix(int ** substitution_matrix) {
	int size = aminoacid_alphabetSize+1;


	cout << endl;
	cout << "\t" ;
	for (int j=0; j<size; j++) {
		cout << aminoacid_num2ascii[j] << "\t";
	}
	cout << endl;
	for (int i=0; i<size; i++) {
		cout << aminoacid_num2ascii[i] << "\t";
		for (int j=0; j<size; j++) {
				cout << substitution_matrix[i][j] << "\t";
		}
		cout << endl;
	}

}


void Class_Substitution_Matrix::checkMatrix(int ** substitution_matrix) {
	int size = aminoacid_alphabetSize+1;


	for (int i=0; i<size; i++) {
		for (int j=0; j<size; j++) {
				if (substitution_matrix[i][j] == INT_MAX) {
						cerr << "matrix field not filled !?" << endl;
						exit(1);
				}
		}
	}
}

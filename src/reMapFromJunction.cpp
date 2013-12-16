#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <set>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <cstdlib>

#define _USE_MATH_DEFINES
using namespace std;

set<string> total_chr;

map<string, set<int> > firstPoint, secondPoint, lockFirst;

map<string, map<int,set<int> > > Nbor;

map<string, map<int, set<int> > > forJunc, revJunc;

int Len;//read length

string single_ref;//ref

void search_Nbor(vector<int>&, int, string);
void print_new_ref(string);
void findNearNbor();
pair<int, int>junc;
map<string, map< pair<int, int>, double > > juncScore;



int run_from_ref(ifstream& input_ref){
	string fa_line,chr,chr_temp;
	int flag = 0;
	int F = 1;
	while(!input_ref.eof()){
		getline(input_ref, fa_line);
		if(fa_line.length() < 1)
			continue;
		if(fa_line.substr(0,1) == ">"){
			istringstream ss(fa_line);
			string fa;
			ss >> fa;
			chr = fa.substr(1);
			if(flag == 0){
				flag = 1;
				chr_temp = chr;
			}
			else{
				if(F == 1){
					transform(single_ref.begin(),single_ref.end(),single_ref.begin(),::toupper);
					print_new_ref(chr_temp);
					total_chr.erase(chr_temp);
				}
				if(chr_temp.size() == 0)
					return 0;
				chr_temp = chr;
				single_ref = "";
			}

			if(total_chr.find(chr) == total_chr.end())
				F = 0;
			else
				F = 1;
			continue;

		}
		if(F == 1)
			single_ref += fa_line;
	}
	transform(single_ref.begin(),single_ref.end(),single_ref.begin(),::toupper);
	print_new_ref(chr_temp);
	return 0;
}

void getJuncs(ifstream& input_junc){
	string junc_line, chr , temp_b, temp_e, temp_s,t;
	int b, e;
	double score;
	while(!input_junc.eof()){
		getline(input_junc, junc_line);
		istringstream ss(junc_line);
		ss >> chr >> temp_b >> temp_e >> t >> temp_s;
		total_chr.insert(chr);
		b = atoi(temp_b.c_str());
		e = atoi(temp_e.c_str());
		score = atof(temp_s.c_str());
		if(score >= 0.5){
		firstPoint[chr].insert(b);
		secondPoint[chr].insert(e);
		}
		forJunc[chr][b].insert(e);
		revJunc[chr][e].insert(b);
		junc.first = b;
		junc.second = e;
		juncScore[chr][junc] = score;

	}
	findNearNbor();
}

void findNearNbor(){
	for(map<string, set<int> >::iterator it = secondPoint.begin() ; it != secondPoint.end(); it++){
		for(set<int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			set<int>::iterator it3 = firstPoint[it->first].upper_bound(*it2);
			while(it3 != firstPoint[it->first].end() && *it3 - *it2 < Len){
				Nbor[it->first][*it2].insert(*it3);
				//lockFirst[it->first].insert(*it3);
				it3++;
			}
		}
	}

}

void print(vector<int>& posList, string chr){
	double score = 1;
	for(int i = 0; i <= (posList.size()-2)/2; i++ ){
		junc.first = posList[i*2];
		junc.second = posList[i*2+1];
		score *= juncScore[chr][junc];
	}
	if(score < 0.9 && (posList.size()-2)/2 > 0){
		return;
	}

	int b = max(posList[0]-Len-1,0);
	string temp = single_ref.substr(b, posList[0]-b-1);
	stringstream ss;
	ss << ":" << b+1 << ":" << posList[0]; 
	for(int i = 1; i <= (posList.size()-2)/2; i++ ){
		temp += single_ref.substr(posList[i*2-1]-1, posList[i*2] - posList[i*2-1]);
		ss << ":"<<posList[i*2-1] <<":"<<posList[i*2];
	}
	int e = min(int(single_ref.length()), posList[posList.size()-1]+Len-1);
	ss << ":"<< posList[posList.size()-1]  <<":"<<e << ":" << score;
	temp += single_ref.substr(posList[posList.size()-1]-1, e - posList[posList.size()-1]-1);
	string head = ">" + chr + ss.str();
	cout << head << endl;
	cout << temp << endl;
}

void search_Nbor(vector<int>& posList, int e, string chr, int temp_len){
	if(Nbor[chr].find(e) != Nbor[chr].end() && temp_len < Len){
		for(set<int>::iterator it = Nbor[chr][e].begin(); it != Nbor[chr][e].end(); it++){
			int b = *it;
			int temp_len_new = temp_len + *it - e;
			posList.push_back(b);
			for(set<int>::iterator it3 = forJunc[chr][b].begin(); it3 != forJunc[chr][b].end(); it3++){
				posList.push_back(*it3);
				search_Nbor(posList, *it3, chr, temp_len_new);
				posList.pop_back();
			}
			posList.pop_back();
		}
	}
	else
		print(posList, chr);

}






void print_new_ref(string chr){
	for(map<int, set<int> >::iterator it = forJunc[chr].begin(); it != forJunc[chr].end(); it++){
		for(set<int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			vector<int> posList; 
			posList.push_back(it->first);
			posList.push_back(*it2);
			print(posList, chr);
			if(Nbor[chr].find(*it2) != Nbor[chr].end()){// && lockFirst[chr].find(it->first) == lockFirst[chr].end()){
				search_Nbor(posList, *it2, chr, 0);
			}


		}
	}
}


int main(int argc, char *argv[]){//junction file || reference || Len of reads || candidate reads(fastq)
	ifstream input_junc(argv[1]);
	Len = atoi(argv[3]);
	getJuncs(input_junc);
	ifstream input_ref(argv[2]);
	run_from_ref(input_ref);
}


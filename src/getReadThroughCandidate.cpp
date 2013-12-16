#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <set>
#include <algorithm>
#include <sstream>
using namespace std;

int pairNum;
struct interval{
     long start, end; // [start, end]
     interval(long  s, long  e) : start(s), end(e) {}
};

bool operator<(interval const& a, interval const& b){
     return a.start < b.start;
}

typedef map<interval,set<string> > intervalMap;
typedef map<string, intervalMap> CHRintervalMap;


CHRintervalMap ExonList,IntronList,Gene; // chr pair name
CHRintervalMap::iterator ItExonList,ItGene;
intervalMap Exon;
intervalMap::iterator ItExon;
//interval::iterator Itinterval;
set<string>::iterator ItGeneName,ItGeneName1,ItGeneName2;
pair<string,string> genePair;
map< pair<string,string>, set<string> > readTh;

map<string, pair<long,long> > Genelist;
map<string, pair<long,long> >::iterator ItGenelist;

map<string, string> GeneChr;
map<string, string> GeneOri;

//				cout <<"\t"<<chr2<<"\t"<<begin2;

interval const*find(intervalMap& s, long  point)
{
     intervalMap::iterator Itinterval = s.lower_bound(interval(point,point));
     if(Itinterval == s.end() || point < Itinterval->first.start) {
         if(Itinterval == s.begin())
             return NULL;
         --Itinterval;
         return point <= Itinterval->first.end && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;
     }
     return &(Itinterval->first);

}
/*
void check(interval_set & s,  point)
{
     if(interval *i = find(s, point))
         printf("%f belongs to [%f,%f]\n", point, i->start, i->end);
     else
         printf("%f does not belong to any interval\n", point);

} 

*/

void merge_interval(intervalMap* map, long  b, long  e, string Name){
	if((*map).find(interval(b,e)) != (*map).end()){
		if((*map)[interval(b,e)].find(Name) !=(*map)[interval(b,e)].end())
			return;
		else{
			(*map)[interval(b,e)].insert(Name);
			return;
		}
	}
	else{
//		interval *i = find(map, b);
//		int end = e > i->first->end? e:i->first->end;
		(*map)[interval(b,e)].insert(Name);
	}
}
		

interval const* searchName(intervalMap* map, long  b, long  e){
	intervalMap::iterator Itinterval = (*map).lower_bound(interval(b,e));
	if(Itinterval == (*map).end() || b < Itinterval->first.start) {
                if(Itinterval == (*map).begin())
                return NULL;
         	--Itinterval;
                return b <= Itinterval->first.end && b >= Itinterval->first.start ? &(Itinterval->first) : NULL;
        }
        return &(Itinterval->first);

}


interval const* rt(intervalMap* map, long  b, long b1){
        intervalMap::iterator Itinterval = (*map).lower_bound(interval(b,b));
        if(Itinterval == (*map).end() || b < Itinterval->first.start) {
                if(Itinterval == (*map).begin())
                return NULL;
                --Itinterval;
                return Itinterval->first.start > b1 && b >= Itinterval->first.start ? &(Itinterval->first) : NULL;
        }
        return Itinterval->first.start > b1 && b >= Itinterval->first.start ? &(Itinterval->first) : NULL;

}



void getGeneList(ifstream& input_gene){
	string gene_line;
	while(!input_gene.eof()){
	        getline(input_gene, gene_line);
		istringstream ss(gene_line);
        	string S,chr,begin,end,name,temp;
        	ss >> S >> chr >> begin >> end >> name;
		long  B = atoi(begin.c_str());
		long  E = atoi(end.c_str());
		if(S == "@"){
			merge_interval(&ExonList[chr],B,E,name);
		
		
			if( Genelist.find(name) == Genelist.end() ){
				 Genelist[name]=pair<long,long>(B,E);
			}
			else{
				Genelist[name].first = Genelist[name].first<B?Genelist[name].first:B;
				Genelist[name].second = Genelist[name].second>E?Genelist[name].second:E;
			}
		}
		if(S == "-"){
                merge_interval(&IntronList[chr],B,E,name);
                }
		GeneChr[name] = chr;
	}
}


void getGeneInterval(){
	for(ItGenelist=Genelist.begin();ItGenelist!=Genelist.end();ItGenelist++){
		merge_interval(&Gene[GeneChr[ItGenelist->first]],ItGenelist->second.first,ItGenelist->second.second,ItGenelist->first);
	}
}


void checkPairReads(ifstream& input_pair){
	string pair_line;
	int flag = 0;
        while(!input_pair.eof()){
                getline(input_pair, pair_line);
                istringstream ss(pair_line);
                string chr1,chr2,begin1,begin2,end1,end2,name, ori1,ori2;
                ss >> name >> chr1 >> begin1 >> end1 >> ori1 >> chr2 >> begin2 >> end2 >> ori2;
//		cout << name<<"\t"<<begin1<<endl;
		if(chr1!=chr2)continue;

		interval const*F1 = searchName(&Gene[chr1],atoi(begin1.c_str()),atoi(end1.c_str()));
		interval const*F2 = searchName(&Gene[chr2],atoi(begin2.c_str()),atoi(end2.c_str()));
	//	if(F1==NULL)F1 = searchName(&IntronList[chr1],atoi(begin1.c_str()),atoi(end1.c_str()));
	//	if(F2==NULL)F2 = searchName(&IntronList[chr2],atoi(begin2.c_str()),atoi(end2.c_str()));
		flag = 0;
		if(F1!=NULL&&F2!=NULL){
//			cout << name << "\t"<<(*F1).start << "\t"<<(*F1).end<<"\t"<<begin1<<"\t"<<(*F2).start<<"\t"<<(*F2).end<<"\t"<<begin2<<endl;

			for(ItGeneName1 = (Gene[chr1])[*F1].begin(); ItGeneName1 != (Gene[chr1])[*F1].end(); ItGeneName1++){
				for(ItGeneName2 = (Gene[chr1])[*F2].begin();ItGeneName2 != (Gene[chr1])[*F2].end(); ItGeneName2++){
					if( *ItGeneName1 == *ItGeneName2){
						flag = 1;
						break;
					}
				}
				if(flag == 1)
					break;
			}
			if(flag == 0){
//				cout << name << "\t"<< chr1;
				for(ItGeneName1 = (Gene[chr1])[*F1].begin(); ItGeneName1 != (Gene[chr1])[*F1].end(); ItGeneName1++){
					
	                                for(ItGeneName2 = (Gene[chr1])[*F2].begin();ItGeneName2 != (Gene[chr1])[*F2].end(); ItGeneName2++){
						if(Genelist[*ItGeneName2].first>Genelist[*ItGeneName1].first){
							genePair.first = *ItGeneName1;
							genePair.second = *ItGeneName2;
						}
						else{
							genePair.first = *ItGeneName2;
                                                        genePair.second = *ItGeneName1;
						}
//						 cout << "\t"<< *ItGeneName1<<"\t"<<*ItGeneName2<<endl;
					}
					readTh[genePair].insert(name);
				}
			}

		}
	//	flag = 0;
	}
	ofstream fout("Namelist");
	for(map< pair<string,string>, set<string> >::iterator it = readTh.begin();it!=readTh.end();it++){
		if(int(it->second.size()) < int(pairNum)){continue;}
		if(GeneOri[it->first.first] != GeneOri[it->first.second]){continue;}
		long f1,f2,e1,e2;
		f1 = Genelist[it->first.first].first;
		e1 = Genelist[it->first.first].second;
		f2 = Genelist[it->first.second].first;
		e2 = Genelist[it->first.second].second;
		if( ( f1 <= e2 && f1 >= f2 ) || ( e1 <= e2 && e1 >= f2 ) || (f2 <= e1 && f2 >= f1) ){continue;}
//		cout << it->first.first <<"\t"<<it->first.second<<"\t"<<it->second.size()<<"\t";
//		for(set<string>::iterator itset = it->second.begin(); itset != it->second.end();itset++){
//			cout << "\t"<<*itset;
//		}
		
//		cout << GeneChr[it->first.first]<<"\t"<<Genelist[it->first.first].first<<"\t"<<Genelist[it->first.first].second<<"\t"<<Genelist[it->first.second].first<<"\t"<<Genelist[it->first.second].second<<endl;
	//	if(f2-e1>40000){continue;}
		interval const*F1 = rt(&Gene[GeneChr[it->first.first]],f2-10,e1);
		if(F1 == NULL || f2-e1<40000){
			fout << it->first.first << "\t"<<it->first.second << endl;
			printf("%s:%ld-%ld\t%s:%ld-%ld\t[+ +]\t(%d)\n",GeneChr[it->first.first].c_str(),Genelist[it->first.first].first,Genelist[it->first.first].second,GeneChr[it->first.first].c_str(),Genelist[it->first.second].first,Genelist[it->first.second].second,int(it->second.size()));
		}
	}
}
			
/*
void getCandidate(ifstream& input_region){
	string region_line;
	while(!input_region.eof()){
                getline(input_pair, pair_line);
                istringstream ss(pair_line);
                string chr1,chr2,begin1,begin2,end1,end2,name, ori1,ori2;
                ss >> name >> chr1 >> begin1 >> end1 >> ori1 >> chr2 >> begin2 >> end2 >> ori2;



}
*/
void Print(){
	 for(ItExonList = ExonList.begin(); ItExonList != ExonList.end(); ItExonList++){
//		cout << ItExonList->first<<"\t";
		for(ItExon = ItExonList->second.begin(); ItExon != ItExonList->second.end(); ItExon++){
			cout << ItExonList->first<<"\t";
			cout << ItExon->first.start <<"\t"<< ItExon->first.end;
			for(ItGeneName = ItExon->second.begin(); ItGeneName != ItExon->second.end();  ItGeneName++)
				cout << "\t"<< *ItGeneName;
			cout<<endl;
		}
	}
}

void getKnownGene(ifstream& input_known){
	string known_line;
	getline(input_known, known_line);
        while(!input_known.eof()){
                getline(input_known, known_line);
		istringstream ss(known_line);
                string id,id2,chr,temp1,temp2,name, ori;
		int t1,t2,t3,t4,t5;
                ss >> id >> chr >> ori  >> t1>>t2>>t3>>t4>>t5>>temp1 >> temp2 >> id2 >> name;
	//	sscanf(known_line.c_str(),"%*s %*s %s %*d %*d %*d %*d %*d %*s %*s %s",ori,Name);
	//	string name = Name;
	//	cout << name <<"\t"<<ori<<endl;
		GeneOri[name] = ori;
	}
}


int main(int argc, char *argv[]){
	ifstream input_gene(argv[argc-1]);
	ifstream input_pair(argv[argc-2]);
	ifstream input_known(argv[argc-3]);
	pairNum = atoi(argv[argc-4]);
	getKnownGene(input_known);
//	ifstream input_region(argv[argc-3]);
	getGeneList(input_gene);
	getGeneInterval();
//	Print();
	checkPairReads(input_pair);
	
	
	return(0);
}
	
		

		
	

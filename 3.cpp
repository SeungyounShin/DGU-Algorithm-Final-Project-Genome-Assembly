#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <map>
#include <algorithm>
using namespace std;
//2018112005 컴퓨터공학과 신승윤

void init(int length){
  //create RNA seq
  ofstream file;
  file.open("RNA_seq.txt");
  srand(time(NULL));
  for(int i = 0; i<length; i++){
      int randNum = rand() % 4;
      if(randNum == 0) file << 'A';
      else if(randNum == 1) file << 'C';
      else if(randNum == 2) file << 'G';
      else file << 'T';
  }
  file.close();
}

string import(){
  string line;
  ifstream file("RNA_revise.txt");
  getline (file,line);
  file.close();
  return line;
}

void reviseSeq(){
  string line;
  srand(time(NULL));
  ifstream file("RNA_seq.txt");
  getline (file,line);
  for(int i=0; i<line.length()-1; i++){
    if(line[i] == line[i+1]){
      int randNum = rand() % 4;
      char c;
      if(randNum == 0) c='A';
      else if(randNum == 1) c ='C';
      else if(randNum == 2) c='G';
      else c='T';
      line[i] = c;
    }
  }
  file.close();
  ofstream file2;
  file2.open("RNA_revise.txt");
  file2 << line;
  file2.close();
  return ;
}

vector<string> Sequencer(int n, int k, string s){
  vector <string> rtn;
  vector <int> rands;

  random_device rd1,rd2,rd3,rd4;
  mt19937 mt(rd1());
  mt19937 mt2(rd2());
  mt19937 mt3(rd3());
  mt19937 mt4(rd4());
  uniform_int_distribution<int> dist(0, s.length()-k-2);
  uniform_int_distribution<int> dist2(0, 2);
  uniform_int_distribution<int> dist3(0, k-1);
  uniform_int_distribution<int> dist4(0, 4);

  //make short read
  bool check[s.length()];
  char myseq[s.length()];
  for(int i=0; i<s.length(); i++) check[i] = false;
  for(int i=0; i<s.length(); i++) myseq[i] = 'x';
  for(int i=0; i<n; i++){
    int randNum = dist(mt);
    int cover = 0;
    for(int j=0; j<k ;j++) if(check[randNum + j]) cover+=1;
    if(cover >= k-10){i--; continue;}//for uniform distributed reads
    rands.push_back(randNum);
    string sub = s.substr(randNum, k);
    for(int j=0; j<k; j++){
      myseq[randNum+j] = sub[j];
      check[randNum+j] = true;
    }
    //for(int j=0; j<s.length() ;j++) cout << check[j] ;
    //cout << endl;
    rtn.push_back(sub);
  }

  cout << "[shortReads maked] :: "<< rtn.size() <<endl;
  int matched = 0;
  int is_uniform = 0;
  for(int i=0; i<s.length(); i++){
    if(s[i] == myseq[i]) matched += 1;
    if(check[i]) is_uniform+=1;
  }
  cout << " matched :: " << matched << "/" << s.length() << endl;
  cout << " cover   :: " << is_uniform << "/" << s.length() << endl;
  // save to txt
  ofstream file;
  file.open("shortReads.txt");
  // noise (SNIP) to short reads
  for(auto read : rtn){
    int noiseNum = dist2(mt2);
    for(int j=0; j<noiseNum; j++){
      int randInd = dist3(mt3);
      int randChar = dist4(mt4);
      char c;
      if(randChar == 0) c='A';
      else if(randChar == 1) c='C';
      else if(randChar == 2) c= 'G';
      else c='T';
      read[randInd] = c;
    }
    file << read << endl;
  }
  file.close();

  ofstream file2;
  file2.open("rands_debug.txt");
  for(auto read : rands){
    file2 << read << endl;
  }
  file2.close();

  ofstream file3;
  file2.open("myseq_debug.txt");
  file2 << string(myseq) << endl;
  file2.close();
  return rtn;
}

void BruteForceRecon(string ref, vector<string> &shortReads){
  // reconstruct DNA
  // init
  clock_t start = clock();
  int k = shortReads[0].length();
  int n = shortReads.size();
  char recon[ref.length()];
  cout << "init done!" << endl;
  for(int i=0; i<ref.length(); i++) recon[i]='A';
  // matching
  for(int readInd =0; readInd<shortReads.size(); readInd++){
    if(readInd % 1000 == 0) cout << readInd << " progress!" << endl;
    string read = shortReads[readInd];
    for(int j=0; j<ref.length()-k-1; j++){
      //sliding window
      int diff = 0;
      // diff calculation
      for(int ind=0; ind<k; ind++){
        if(diff > 2) break;
        if(read[ind] != ref[j+ind]){
          diff += 1;
        }
      }
      // match
      if(diff <= 2){
        for(int i=0; i<k; i++) recon[j+i] = read[i];
      }
    }
  }
  // matching accuracy
  int matched = 0;
  for(int i=0; i<ref.length(); i++){
    if(ref[i] == recon[i]) matched += 1;
  }
  clock_t end = clock();
  double duration = (double)(end - start) / CLOCKS_PER_SEC;
  cout << " Time spend        :: " << duration << endl;
  cout << " matching accuracy :: " << matched << "/" << ref.length() << endl;
  return ;
}

string ShortestCommonSuperstring(vector<string> &shortReads){
  string shortest_sup;
  int n = shortReads.size();
  vector<int> perm(n);
  vector<int> start(n);
  for(int i=0; i<n; i++) perm[i]=i;
  for(int i=0; i<n; i++) start[i]=i;
  while (next_permutation(perm.begin(), perm.end())) {
    for (auto& i : perm)
        cout << i << "  ";
        cout << '\n';
  }
}

void denovo(vector<string> &shortReads){
  cout << "denovo" << endl;
  map <string, vector<string> > G;
  int k = shortReads[0].length();
  for(auto read : shortReads){
    // read :: "ATC...TGC"
    string left = read.substr(0,k-1);
    string right = read.substr(1,k-1);

    G[left].push_back(right);

  }

  // Graph Print
  for(auto read : shortReads){
    string left = read.substr(0,k-1);
    string right = read.substr(1,k-1);

    cout << left << " -> ";
    for(int i=0; i<G[left].size(); i++){
      cout << G[left][i] << " ";
    }
    cout << endl;
  }


  return;
}

int main(){
  //init(1000000); // random Sequence generate
  //reviseSeq(); // duplicate removal

  string ref = import();
  //vector<string> shortReads = Sequencer(40000, 30, ref);
  //string arr[] = {"a_lon", "_long", "long_", "ong_l","ng_lo","g_lon",
  //                "_long","long_","ong_t","ng_ti","g_tim","_time"};
  string arr[] = {"ACGGATGAGC", "GAGCGGA", "GAGCGAG"};
  // Initialize vector with a string array
  vector<string> shortReads(arr, arr + sizeof(arr)/sizeof(string));
  //BruteForceRecon(ref, shortReads);
  cout << ShortestCommonSuperstring(shortReads) <<endl;

  return 0;
}

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <map>
#include <algorithm>
#include <tuple>
using namespace std;
//2018112005 컴퓨터공학과 신승윤

string import(string path){
  string line;
  ifstream file(path);
  getline (file,line);
  //cout << line << endl;
  file.close();
  return line;
}

void ShortRead(int L, int M, string s) { //shortread 생성
    ofstream fout;
    ifstream fin;
    string shortread;
    int idx;
    srand(time(NULL));
    fout.open("shortread_" + to_string(L) + "_" + to_string(M) + ".txt"); //해당 파일 열거나 없으면 생성
    for (int i = 0; i < M; i++) {
        idx = rand() % (s.length() - L); //랜덤 위치 얻기
        shortread = s.substr(idx, L);
        fout << shortread << endl;
        shortread = "";
    }
    fin.close();
    fout.close(); //파일닫기
}

vector<string> importReads(string path, int L, int M){
  vector<string> ret;
  ifstream file(path);

  int count = 0;
  if(file.is_open()){
    while(!file.eof()){
      string tmp;
      if(ret.size() == M) return ret;
      getline(file,tmp);
      ret.push_back(tmp.substr(0,L));
      count++;
    }
  }
  file.close();    //파일 닫기
  return ret;
}

int overlap(string s1, string s2,int min_length){
  int start = 0;
  cout << " overlap func :: "<< s1 << " " << s2 << endl;
  start = s1.find(s2.substr(0,min_length));
  cout << start << " "<< s1 << " "<<s2<<endl;
  if(start == -1){
    if(s1.length() >= s2.length()){
      for(int i=1;i<s2.length();i++)
        if(s1.substr(s1.length()-s2.length()+i, s2.length()-i)  == s2.substr(0,s2.length()-i)){
          return s1.length()-s2.length()+i;
        }
    }
    else{
      for(int i=1;i<s1.length();i++)
        if(s1.substr(i,s1.length()-i)  == s2.substr(0,s1.length()-i)){
          return i;
        }
      //is s1 is subset of s2?
    }
    return start;
  }
  return start;
}

int overlap2(string s1, string s2,int min_length){
  int start = -1;
  //case 1 :: full overlap
  // ------
  // ------
  start = s1.find(s2);
  if(start > 0 && start < s1.length() && s1.length() > s2.length()){
    return start;
  }
  //case 2 :: ------
  //              --------
  for(int i=1; i<s1.length(); i++){
      string right = s1.substr(i,s1.length()-i);
      string left = s2.substr(0,s1.length()-i);
      //cout << "CASE 2 " << endl;
      //cout << s1 << endl;
      //for(int j=0; j<i;j++) cout << " "; cout << s2 << endl;
      if(right == left){
        return i;
      }
  }

  return start;
}

tuple <int,int,int> pick_maximal_overlap(vector <string> reads, int k, bool prnt){
  bool print = prnt;
  vector <int> out;
  int best_olen = 0;
  int reada, readb;

  if(print) cout << " pick maximal overlap :: reads " << endl;
  if(print) for(auto s : reads)cout << " " << s << endl;
  if(print) cout << "=====" << endl;

  map <string, vector<string> > index;
  int cnt = 0;
  for(auto read : reads){
    for(int i=0; i<read.length()-k+1; i++){
      index[read.substr(i,k)].push_back(read);
      //if(print) cout << read.substr(i,k) << " ";
    }
    cnt ++;
    //if(print) cout << endl;
  }
  //cout<<"flag"<<endl;

  int i,j;i=0;j=0;
  int olenInd = 0;
  int count = 0;
  string tmpa; string tmpb;
  for(auto r : reads){
    for(auto o : index[r.substr(r.length()-k,k)] ){
      count += 1;
      //cout<< r.length() << " " << o.length() << endl;
      if(print) cout << " r :: "<< r << " o ::" << o << " ||"<< endl;
      //if(print)cout << " sub min : " << r.substr(r.length()-k,k) << endl;
      //if(print)for(auto sub : index[r.substr(r.length()-k,k)]) cout << "    " << sub<<endl;
      if(r != o){
        bool flag = false;
        int olen = overlap2(r,o,k);
        int overlap = r.length() - olen;
        int min_len = min(o.length(), r.length());
        if(olen == -1) overlap = 0;
        //if(print) cout << " + overlap " << overlap<< " > " << best_olen << " "<< olen << " "<< r << " "<<o<< " " << endl;

        if(overlap >= min_len){
          tmpa = r; tmpb = o;

          //cout << r << " "<< o << endl;
          best_olen = overlap;
          olenInd = olen;
          break;
        }
        if(overlap > best_olen){
          tmpa = r; tmpb = o;

          //if(print) cout << "best overlap " << overlap<< " "<< r << " "<<o<< " " << endl;
          best_olen = overlap;
          olenInd = olen;
        }
      }
    }
  }
  vector <string>::iterator itera = find(reads.begin(), reads.end(), tmpa);
  vector <string>::iterator iterb = find(reads.begin(), reads.end(), tmpb);
  reada = distance(reads.begin(), itera);
  readb = distance(reads.begin(), iterb);
  //cout << " rept :: " << count << endl;
  return make_tuple(reada, readb, olenInd);
}

vector <string> greedySCS(vector <string> shortReads, int k){
  bool print =false;
  clock_t start = clock();
  int M = shortReads.size();
  int L = shortReads[0].length();
  tuple<int,int,int> out = pick_maximal_overlap(shortReads, k, false);
  int reada = get<0>(out); int readb = get<1>(out); int olen =get<2>(out);
  int count = 0;
  cout << shortReads[reada] << " " << shortReads[readb] << " "<<olen << endl;
  while(olen>-1 && shortReads.size() > 1){
    // adaptive scheduling
    if(count < M/10){
      k = L-1;
    }else if(count >= M/10 && count <M/2){
      k = L-10;
    }else{
      k = L/2;
    }
    cout << "progress :: "<<count << " "<< shortReads.size() << endl;
    if(print) cout << " [shortReads-before-erase] " << shortReads.size() <<endl;
    if(print) for(auto t : shortReads) cout << " " << t << endl;
    if(print) cout << "greedySCS - choose "<<  reada << " "<< readb << " "<< olen << endl;
    if(reada == readb) break;
    string read_a = shortReads[reada];
    string read_b = shortReads[readb];
    if(print) cout << "     "<< read_a << " " << read_b << endl;
    if(reada > readb){
      shortReads.erase(shortReads.begin()+reada);
      shortReads.erase(shortReads.begin()+readb);
    }else{
      shortReads.erase(shortReads.begin()+readb);
      shortReads.erase(shortReads.begin()+reada);
    }
    if(print) cout << " [shortReads-after-erase] " << shortReads.size() <<endl;
    if(print) for(auto t : shortReads) cout << " " << t << endl;
    if(read_a.find(read_b) > 0 &&read_a.find(read_b) <= read_a.length() && (read_a.length() > read_b.length())){
      if(print) cout << " push[1] -> " << read_a << " "<<  read_a.find(read_b)<<endl;
      if(print) cout << "       -> " << read_a << " "<<  read_b<< endl;
      shortReads.push_back(read_a);
    }else if(read_b.find(read_a) > 0 &&read_b.find(read_a) <= read_b.length() && (read_b.length() > read_a.length())){
      if(print) cout << " push[2] -> " << read_b << " "<<  read_b.find(read_a)<< endl;
      if(print) cout << "       -> " << read_a << " "<<  read_b<< endl;
      shortReads.push_back(read_b);
    }else{
      if(print) cout << " push[3] -> " << read_a.substr(0, olen) << " + " << read_b << endl;
      shortReads.push_back(read_a.substr(0, olen) + read_b);
    }
    if(print) cout << " [shortReads-after-insert] " << shortReads.size() <<endl;
    if(print) for(auto t : shortReads) cout << " " << t << endl;
    tuple<int,int,int> tmp = pick_maximal_overlap(shortReads, k, false);
    reada = get<0>(tmp); readb = get<1>(tmp); olen =get<2>(tmp);
    count += 1;
  }
  string ret = "";
  clock_t end = clock();
  cout << "Time spend :: "<< (end-start)/1000000 <<endl;
  return shortReads;
}

int main(){
  int M, L,N;
  N = 500000;
  L = 100;
  M = 110000;
  string DIR_PATH = "./denovo/" + to_string(N) + "/" + to_string(L) + "_"+ to_string(M) +"/";
  string ref = import(DIR_PATH + "mydna.txt");
  N = ref.length();
  vector<string> shortReads = importReads(DIR_PATH + "shortread.txt",L, M);

  cout << " # of ShortReads :: "<< shortReads.size() << endl;
  //string ref = "AGATAACTGGGCCCTACGGCAATATACGTACAATTTAGGA";
  //string arr[] = {"ACGGATGAGC", "GAGCGGA", "GAGCGAG"};
  //string arr[] = {"AGATAACTG", "TAACTGGGC", "TAACTGGGCCC","AACTGGGC","TGGGCCCTACG","CCCTACGGCAA","CTACGGCAATATAC","CGGCAATATACGT","ATATACGTACAATTT","CGTACAATTTAGGA"}; // Initialize vector with a string array
  //vector<string> shortReads(arr, arr + sizeof(arr)/sizeof(string));

  vector<string> recon = greedySCS(shortReads, L-1);
  cout << " + recon size :: "<<recon.size() << endl;

  ofstream file("ref.txt");
  ofstream file2("recon_vec.txt");
  ofstream file3("shortreads.txt");
  file << ref;
  for(auto r : recon){
      file2 << r << endl;
  }
  for(auto s : shortReads){
      file3 << s << endl;
  }
  file.close();
  file2.close();
  file3.close();

  if(recon.size() == 1){
    int matched = 0;
    for(int i=0; i<N; i++){
      if(ref[i] == recon[0][i]) matched++;
    }
    cout << matched << "/" << N << endl;
  }

  return 0;
}

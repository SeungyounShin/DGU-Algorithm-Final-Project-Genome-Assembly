//2018112017 이다연
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
using namespace std;
#define D 3 //미스매치 수




vector<vector<int>> filltable(int k, int n, vector<vector<int>> table) {//테이블생성함수
    
    ifstream fin;
    string reference;
    int sumofascii = 0;
    fin.open("Reference_DNA_sequence.txt");
    getline(fin, reference);
    for (int i = 0; i < reference.size(); i++) {
        for (int j = 0; j < k / (D + 1); j++) {
            sumofascii += int(reference[i + j]); //아스키값 총합을 인덱스로
        }
        table[sumofascii].push_back(i); //테이블에 넣음
        //cout << sumofascii << endl;
        sumofascii = 0;
    }
    return table;
}

void perfectMatch(int k, int n, vector<vector<int>> table) { //퍼펙트매칭함수
    ifstream fin1, fin2;
    ofstream fout; //ofstream 형식 변수 선언
    string shortread, reference;
    int tmp = 0;
    int cnt = 0;
    int isExit = 0;
    fin1.open("ShortRead" + to_string(k) + "_" + to_string(n) + ".txt"); //해당 파일 열거나 없으면 생성
    fin2.open("Reference_DNA_sequence.txt");
    fout.open("Reconstruct_DNA_sequence_" + to_string(k) + "_" + to_string(n) + ".txt");
    string reconstruct;
    getline(fin2, reference);
    reconstruct = reference;
    
    for (int i = 0; i < n; i++) { //shortread 개수만큼
        getline(fin1, shortread);
        isExit = 0;
        for (int j = 0; j < (D + 1); j++) { //partition 만큼 반복
            tmp = 0;
            if (isExit == 1) { //퍼펙트매치 이미 발견했을경우
                break;
            }
            for (int l = 0; l < k / (D + 1); l++) { //partition 길이만큼 반복
                //cout << shortread[j * k / (MISSMATCH + 1) + l];
                tmp += int(shortread[j * k / (D + 1) + l]);
            }
            for (int m = 0; m < table[tmp].size(); m++) { //테이블 index 개수만큼
                for (int n = table[tmp][m]; n < table[tmp][m] + k / (D + 1); n++) { //index와 shortread 비교
                    if (shortread[cnt] != reference[n]) {
                        cnt = 0;
                        break;
                    }
                    else {
                        cnt++;
                    }
                }
                if (cnt == k / (D + 1)) { //만약 partition이 perfect match면
                    // cout << shortread << endl;
                    for (int o = 0; o < shortread.size(); o++) { // 복원될 string에 넣음
                        reconstruct[table[tmp][m] + o] = shortread[o];
                        //cout << shortread[o];
                    }
                    
                    cnt = 0;
                    isExit = 1;
                    break;
                }
            }
        }
    }
    
    fout << reconstruct;
    fin1.close();
    fin2.close();
    fout.close(); //파일닫기
}

void MyDNA(string s, char Nucleic[], int k, int n) { //k 마다 1개씩 다르도록 sequence 생성
    srand(time(NULL));
    ofstream fout;
    int change;
    int idx;
    string str;
    random_device rd;  //비결정적 생성기
    mt19937 gen(rd()); //메르센 트위스터 시드 설정
    uniform_int_distribution<> dis(0, s.length() - 1); //0부터 reference 길이만큼으로 분포 설정
    fout.open("My_DNA_sequence_" + to_string(k) + "_" + to_string(n) + ".txt"); //해당 파일 열거나 없으면 생성
    for (int i = 0; i < (s.length() * 0.01); i++) {
        change = dis(gen); //바꿀 자릿수
        idx = rand() % 4; //바꿀 문자 인덱스
        if (s[change] != Nucleic[idx]) { //다르면 문자 바꿔주기
            s[change] = Nucleic[idx];
        }
        else {
            while (s[change] == Nucleic[idx]) { //바꿀 문자와 이미 같으면 다른 문자로 바꿔줌
                idx = rand() % 4;
            }
            s[change] = Nucleic[idx];
        }
    }
    fout << s;
    fout.close(); //파일닫기
}
void ShortRead(int k, int n) { //shortread 생성
    ofstream fout;
    ifstream fin;
    string shortread, s;
    int idx;
    srand(time(NULL));
    fin.open("My_DNA_sequence_" + to_string(k) + "_" + to_string(n) + ".txt"); //해당 파일 열거나 없으면 생성
    getline(fin, s); //str 받아오기
    fout.open("ShortRead" + to_string(k) + "_" + to_string(n) + ".txt"); //해당 파일 열거나 없으면 생성
    for (int i = 0; i < n; i++) {
        idx = rand() % (s.length() - k); //랜덤 위치 얻기
        shortread = s.substr(idx, k);
        fout << shortread << endl;
        shortread = "";
    }
    fin.close();
    fout.close(); //파일닫기
}

void Reconstruct(int k, int n) { //brute force 로 복원하는 함수
    ofstream fout; //ofstream 형식 변수 선언
    ifstream fin1, fin2;  //ifstream 형식 변수 선언
    string shortread, reference;
    int cnt; // miss match count 변수
    fin1.open("ShortRead" + to_string(k) + "_" + to_string(n) + ".txt");
    fin2.open("Reference_DNA_sequence.txt");
    fout.open("Reconstruct_DNA_sequence_" + to_string(k) + "_" + to_string(n) + ".txt");
    getline(fin2, reference);
    string reconstruct;
    reconstruct = reference;
    for (int i = 0; i < n; i++) {
        getline(fin1, shortread);
        for (int j = 0; j < n - k; j++) { //brute force 로 하나하나 비교
            cnt = 0;
            for (int l = 0; l < k; l++) {
                if (shortread[l] != reference[j + l]) { //miss match count
                    cnt++;
                }
            }
            if (cnt <= D) { //miss match D 개까지 허용
                for (int l = 0; l < k; l++) {
                    
                    reconstruct[j + l] = shortread[l];
                }
            }
        }
    }
    fout << reconstruct;
    fin1.close();
    fin2.close();
    fout.close(); //파일닫기
}

void Result(int k, int n) { //일치율 출력
    ifstream fin1, fin2;  //ifstream 형식 변수 선언
    string mydna, reconstruct;
    int cnt = 0; //틀린 개수 count
    fin1.open("My_DNA_sequence_" + to_string(k) + "_" + to_string(n) + ".txt"); //해당 파일 열거나 없으면 생성
    fin2.open("Reconstruct_DNA_sequence_" + to_string(k) + "_" + to_string(n) + ".txt");
    getline(fin1, mydna);
    getline(fin2, reconstruct);
    for (int i = 0; i < mydna.length(); i++) {
        if (mydna[i] != reconstruct[i]) {
            cnt++;
        }
    }
    
    printf("일치율: %.2f%%\n", (float(mydna.length()) - float(cnt)) / float(mydna.length()) * 100);
    fin1.close();
    fin2.close();
}

int main() {
    int M, L; // M number of short reads of length L
    clock_t start, end;
    string str; //스트링 변수 선언
    ofstream fout; //ofstream 형식 변수 선언
    ifstream fin;  //ifstream 형식 변수 선언
    cin >> L >> M;
    char Nucleic[4] = { 'A','G','C','T' };
    fout.open("Reference_DNA_sequence.txt"); //해당 파일 열거나 없으면 생성
    random_device rd;  //비결정적 생성기
    mt19937 gen(rd()); //메르센 트위스터 시드 설정
    uniform_int_distribution<> dis(0, 3); //0부터 3으로 분포 설정
    for (int i = 0; i < 1000000; i++) { //숫자만큼 반복하여 랜덤 수 생성한 뒤 문자열에 붙임
        str += Nucleic[dis(gen)];
    }
    fout.write(str.c_str(), str.size()); //문자열 파일에 쓰기
    fout.close(); //파일 닫기
    
    string sequence;  //문자열 저장할 변수 선언
    fin.open("Reference_DNA_sequence.txt"); //시퀀스 파일 열기
    getline(fin, sequence); //str 받아오기
    MyDNA(sequence, Nucleic, 30, 40000); //MyDNA 생성
    //   MyDNA(sequence, Nucleic, 60, 30000);
    ShortRead(30, 40000); //shortread 생성
    //   ShortRead(60, 30000);
    fin.close(); //파일닫기
    
    cout << "Trivial Algorithm" << endl; //브루트포스
    start = clock();
    Reconstruct(L, M); //복원
    end = clock();
    cout << "복원하는 데 걸리는 시간 : " << (double)(end - start) / CLOCKS_PER_SEC << "s" << "\n" << endl;
    Result(L, M);
    
    cout << endl;
    
    vector<vector<int>> table;
    table.resize(1000001);
    cout << "Perfect Matching Algorithm" << endl; //퍼펙트매칭
    start = clock();
    table = filltable(L, M, table); //테이블채우기
    perfectMatch(L, M, table); //매칭하기
    end = clock();
    cout << "복원하는 데 걸리는 시간 : " << (double)(end - start) / CLOCKS_PER_SEC << "s" << "\n" << endl;
    Result(L, M);
    return 0;
}

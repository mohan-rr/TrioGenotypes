#include <iostream>
#include <bitset>

int main(){
    int i;
    int one = 1;
    std::cin >> i;
    std::cout << std::bitset<8*sizeof(i)>(i) << std::endl;
    unsigned char h = i;

    // for (int j = 0; j < sizeof(i); j++){
    for (int k=0; k < 8; k++){
        std::cout << (i & 1);
        i = i >> 1;
        if ( k%2 == 1) std::cout << '\t';
    }
    // std::cout << '\t' ;
    // }
    std::cout << '\n';
    for (int j = 0 ; j < 4; j++){
        std::cout << (h & 3);
        h = h >> 2;
    }
}
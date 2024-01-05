#include <iostream>

int main() {
    std::cout << "Hello, World!" << std::endl;
    
    auto result = (10 <=> 20) > 0;
    if (result) {
        std::cout << "10 > 20" << std::endl;
    } else {
        std::cout << "10 < 20" << std::endl;
    }
}
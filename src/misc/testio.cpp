#include <fstream>
#include <iostream>
 
using namespace std;
 
int main (int argc, char** argv)
{
    // Open a new file for input/output operations discarding any current
    // content in the file (assumes a length of zero on opening)
    fstream myFile("test.txt", ios::in | ios::out | ios::trunc);
 
    // Add the characters "Hello World" to the file
    myFile << "Hello World";
 
    // Seek to 6 characters from the beginning of the file
    myFile.seekg(2, ios::beg);
    myFile.seekg(4, ios::cur);
 
 
    // Read the next 5 characters from the file into a buffer
    char buffer[6];
    myFile.read(buffer, 5);
 
    // End the buffer with a null terminating character
    buffer[5] = 0;
 
    // Output the contents read from the file and close it
    cout << buffer << endl;
 
    myFile.close();
}

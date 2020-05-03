#include <iostream>
#include <iterator> // for std::size
 
bool isVowel(char ch)
{
  switch (ch)
  {
  case 'A':
  case 'a':
  case 'E':
  case 'e':
  case 'I':
  case 'i':
  case 'O':
  case 'o':
  case 'U':
  case 'u':
    return true;
  default:
    return false;
  }
}
 
int main()
{
  char name[]{ "Mollie" };
  int size{std::size(name)};
  std:: cout << size<<'\n';
  int arrayLength{ static_cast<int>(std::size(name)) };
  std::cout<< arrayLength<<'\n';
  int numVowels{ 0 };
 
  for (char* ptr{ name }; ptr < (name + size); ++ptr)
  {
    if (isVowel(*ptr))
    {
      ++numVowels;
    }
  }
 
  std::cout << name << " has " << numVowels << " vowels.\n";
 
  return 0;
}
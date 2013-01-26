/******************************************************************************
File: extractLoci.cpp
Description: The file mloci.out contains multiple loci.out concatenated
             together.  Since multic can only operate on one loci.out at a
             time, we must parse mloci.out into one loci.out.  This program
             does this parsing.
Author: Eric Lunde, 8-23-03
Other: This program is designed to be run from multic.s before each run of the
       alternative hypothesis.  Before the first time this program is run,
       mloci.out is copied to tempmloci.out.  This will allow us to keep the
       information in mloci.out intact while we alter the information in
       tempmloci.out.  After looping in multic.s tempmloci.out should be
       removed.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "verS.h"

extern "C" {

void extractLoci() {
  const int LINE_LENGTH = 256;
  char line[LINE_LENGTH];

  std::ifstream tempmloci("tempmloci.out");
  if(tempmloci.fail()) {
    PROBLEM "The file tempmloci.out cannot be opened for reading.\nextractLoci.cpp key 27\n"
      RECOVER(NULL_ENTRY);
  }

  std::ofstream loci("loci.out");
  if(loci.fail()) {
    PROBLEM "The file loci.out cannot be opened for reading.\nextractLoci.cpp key 35\n"
      RECOVER(NULL_ENTRY);
  }

  // Copy the first line, the narrative line naming the (m)ibd file that
  // created it.
  tempmloci.getline(line, LINE_LENGTH);
  // If tempmloci is already bad, that means that it is empty.
  if(tempmloci.fail()) {
    tempmloci.close();
    loci.close();
    return;
  } 

  loci << line << std::endl;

  // Stop copying when we find '#' or the end of file.
  tempmloci.getline(line, LINE_LENGTH);
  while(tempmloci.good() && line[0] != '#') {
    loci << line << std::endl;
    tempmloci.getline(line, LINE_LENGTH);
  }

  loci.close();

  // Now copy the rest of the file to a temp file.
  std::ofstream deleteThisFile("deleteThisFile");

  while(tempmloci.good()) {
    // line still holds the line of text beginning with '#'
    deleteThisFile << line << std::endl;
    tempmloci.getline(line, LINE_LENGTH);  
  }

  deleteThisFile.close();
  tempmloci.close();

  // Update temploci.out and remove deleteThisFile at the same time
  system("mv deleteThisFile tempmloci.out");
}

}

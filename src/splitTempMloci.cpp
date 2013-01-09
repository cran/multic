#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include "multicString.h"

#include "verS.h"
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif

extern "C" {

/*********************
Title: splitTempMloci
Description: reads a specified mloci.out (typically tempmloci.out) and splits
             it into individual loci.out files (with filenames being that of
             the first line of each loci) in a subdirectory named loci
Input: s_object *fileName - character value specifying the name of the
                            mloci.out file to be read
Output: s_object *ibdNames - character vector whose contents that the ibd or
                             mibd files that genereated a particular loci.out
Side Effects: a directory named 'loci' is created and the loci.out files are
              placed in that directory.  If loci exists already as a
              directory, it and its contents are removed before creating it
              again.
Author: Eric Lunde, 01/12/2005
*********************/
s_object *splitTempMloci(s_object *fileName) {
  S_EVALUATOR

  const int MAX_LINE_LENGTH = 1023;
  char delimiter[MAX_LINE_LENGTH + 1];
  char ibdName[MAX_LINE_LENGTH + 1];
  char line[MAX_LINE_LENGTH + 1];
  char pathName[MAX_LINE_LENGTH + 1];
  char lociDirectory[32] = "loci";
  char *token;

  const char *fileNameValue = CHARACTER_VALUE(fileName);

  // Get the number of loci in this mloci.out file by issuing a
  // 'grep "^#" | wc' command on mloci.out
  char command[MAX_LINE_LENGTH + 1];
  char tempFileName[32] = "splitTempMloci.tmp";
  int lociCount = -1;
  snprintf(command, MAX_LINE_LENGTH, "grep \"^#\" %s | wc > %s",
	   fileNameValue, tempFileName);
  system(command);

  // Read the number of loci from the output file just created by the unix
  // command
  std::ifstream tempFile(tempFileName);
  if(tempFile.fail()) {
    PROBLEM "Error opening %s.\nsplitTempMloci.cpp key 28\n",
      tempFileName RECOVER(NULL_ENTRY);
  }
  tempFile >> lociCount;
  if(lociCount == -1) {
    PROBLEM "Error reading %s\n.slitTempMloci.cpp key 35\n",
      tempFileName RECOVER(NULL_ENTRY);
  }else if(lociCount == 0) {
    return NEW_CHARACTER(lociCount);
  }
  tempFile.close();
  snprintf(command, MAX_LINE_LENGTH, "rm -f %s", tempFileName);
  system(command);

  // Create the s_object to fill with loci names
  s_object *lociNames = NEW_CHARACTER(lociCount);
  int lociNamesIndex = 0;

  // Create the output loci directory
  snprintf(command, MAX_LINE_LENGTH, "test -f %s", lociDirectory);
  if(system(command) == 0) {
    PROBLEM "Cannot create directory '%s' beacuse a file of the same name already exists.\nsplitTempMloci.cpp key 78\n",
      lociDirectory RECOVER(NULL_ENTRY);
  }
  snprintf(command, MAX_LINE_LENGTH, "\\rm -rf %s", lociDirectory);
  system(command);
  snprintf(command, MAX_LINE_LENGTH, "mkdir %s", lociDirectory);
  system(command);
  
  std::ifstream mloci(fileNameValue);
  if(mloci.fail()) {
    PROBLEM "Error opening %s.\nsplitTempMloci.cpp key 44\n",
      fileNameValue RECOVER(NULL_ENTRY);
  }

  // Read the first line of mloci.out
  mloci >> delimiter >> ibdName;
  mloci.ignore(1024, '\n');

  // Create the first loci output file
  snprintf(pathName, MAX_LINE_LENGTH, "%s/%s", lociDirectory, ibdName);
  std::ofstream loci(pathName);
  if(loci.fail()) {
    PROBLEM "Error opening %s.\nsplitTempMloci.cpp key 46\n",
      pathName RECOVER(NULL_ENTRY);
  }

  // Add ibdName to the list of loci names
  multic_SET_STRING_ELT(lociNames, lociNamesIndex, pathName);
  lociNamesIndex++;

  // Begin writing to the first loci file.
  loci << delimiter << " " << ibdName << std::endl;

  // Read mloci until it is empty
  while(mloci.good()) {
    mloci.getline(line, MAX_LINE_LENGTH);

    // If the line is empty, the file is complete
    if(strlen(line) == 0) {
      break;
    }

    // If the line begins with a '#', close this loci file and open a new one
    if(line[0] == '#') {
      loci.close();

      token = strtok(line, " ");
      strcpy(delimiter, token);
      token = strtok(NULL, " ");
      strcpy(ibdName, token);

      // Create the next loci output file
      snprintf(pathName, MAX_LINE_LENGTH, "%s/%s", lociDirectory, ibdName);
      loci.open(pathName);
      if(loci.fail()) {
	PROBLEM "Error opening %s.\nsplitTempMloci.cpp key 65\n",
	  pathName RECOVER(NULL_ENTRY);
      }

      // Add ibdName to the list of loci names
      multic_SET_STRING_ELT(lociNames, lociNamesIndex, pathName);
      lociNamesIndex++;

      // Begin writing to the next loci file.
      loci << delimiter << " " << ibdName << std::endl;

      // Get the next line from mloci
      mloci.getline(line, MAX_LINE_LENGTH);    
    }

    loci << line << std::endl;
  }

  loci.close();
  mloci.close();

  return lociNames;
}

}

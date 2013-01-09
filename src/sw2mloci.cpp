/**
 * sw2mloci is a collection of routines to convert an array of IBD file names
 * and produce one mloci.out.gz file appropriate for multic's immediate use. 
 * The only function of concern is sw2lociFiles.  Call this one from S-PLUS
 * or R.
 *
 * It first converts each IBD file into a series of mibd files.  Then, it
 * combines those mibd files into one mloci.out.gz
 *
 * This could be improved by redesigning the getSwMibdMap* functions.  There
 * is a fair amount of similarity between them and this could EASILY become
 * difficult to maintain.  If a new version of SimWalk comes along (a new
 * version of the IBD files), this redesign would be very useful.
 */
#include <iostream>
#include <fstream>
#include <cctype>
#include <cstring>
#include <map>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>

#include <S.h>
#include "multicString.h"
#include "Rostream.h"
#include "Rstreambuf.h"

using namespace std;
using namespace Rcpp;

// The internal function that can use STL.
void sw2lociFilesInternal(const char **swIbdFileNames,
			  Sint *ibdFileCount,
			  const char **swMapFileName);

// Word/string processing functions.
char *getLastWord(char *line);
bool isEmpty(const char *word);
bool isWordDecimal(const char *line);
bool isWordAlnum(const char *line);
void trimTrailingDecimalZeros(char *decimalNumber);

// IBD- file processing functions.
bool matchesNewPair280(const char *line);
bool matchesNewPair291(const char *line);
bool matchesNewMarker280(const char *line);
bool matchesNewMarker291(const char *line);
bool matchesNewCM280(const char *line);
bool matchesNewCM291(const char *line);

// NameKosambiMap definitions and functions.
typedef pair<const std::string, double> NameKosambiMapEntry;
typedef map<std::string, double> NameKosambiMap;
NameKosambiMap getNameKosambiMap(const char *mapFileName);
void printNameKosambiEntry(NameKosambiMapEntry &ent);
void deleteNameKosambiEntry(NameKosambiMapEntry &ent);

// SwMibdMap definitions and functions.
typedef pair<const double, char *> SwMibdMapEntry;
typedef map<double, char *> SwMibdMap;
SwMibdMap getSwMibdMap(const char *mapFileName,
		       const char *ibdFileName);
SwMibdMap getSwMibdMap280(std::ifstream &ibdFile, const char *mapFileName);
SwMibdMap getSwMibdMap291(std::ifstream &ibdFile, const char *mapFileName);
void printSwMibdEntry(SwMibdMapEntry &ent);
void deleteSwMibdEntry(SwMibdMapEntry &ent);

// CmFileNameMap definitions and functions.
typedef pair<const double, char *> CmFileNameMapEntry;
typedef map<double, char *> CmFileNameMap;
CmFileNameMap getCmFileNameMap(const char *swIbdFileName,
			       const char *swMapFileName,
			       const char *chromosome);
void printCmFileNameEntry(CmFileNameMapEntry &ent);
void deleteCmFileNameEntry(CmFileNameMapEntry &ent);

// CmFileMap definitions and functions.
typedef pair<const double, std::ofstream *> CmFileMapEntry;
typedef map<double, std::ofstream *> CmFileMap;
CmFileMap getCmFileMap(CmFileNameMap &cmFileNameMap, const char *directory);
void printCmFileEntry(CmFileMapEntry &ent);
void deleteCmFileEntry(CmFileMapEntry &ent);
std::ofstream *getMibd(CmFileMap &cm2fileMap, double cM);

// Functions to process the differences in IBD file formats
void ibd2locis280(const char *swIbdFileName, CmFileMap &cmFileMap);
void ibd2locis291(const char *swIbdFileName, CmFileMap &cmFileMap);

// General purpose functions to put similar code in one place.
void lociFiles2mloci(CmFileNameMap &cmFileNameMap, const char *directory);
void addCmsToMap(SwMibdMap &swMibdMap, NameKosambiMap &nameKosambiMap,
		 vector<double> &cMs, const char *previousName,
		 const char *currentName);
void addLastCmToMap(SwMibdMap &swMibdMap, NameKosambiMap &nameKosambiMap,
		    double cM, const char *markerName);
char *int2ascii(int n);
char *double2ascii(double n, int precision);

// Global definitions for string buffer length and text file delimiters.
#define SMALL_BUFFER_LENGTH 32
#define MED_BUFFER_LENGTH 1024
#define TOKEN_DELIMITER " \t\n"
#define DEBUG 1 // PV

/**
 * sw2lociFiles is defined as extern "C" so Splus and R can call it, but I 
 * want to use C++ map structure and cannot use it in a "extern C" function
 * (the Sun CC compiler won't allow it).
 */
extern "C" void sw2lociFiles(const char **swIbdFileNames,
			     Sint *ibdFileCount,
			     const char **swMapFileName) {
  sw2lociFilesInternal(swIbdFileNames, ibdFileCount, swMapFileName);
}

void sw2lociFilesInternal(const char **swIbdFileNames,
			  Sint *ibdFileCount,
			  const char **swMapFileName) {

  Rcout << "Inside sw2lociFilesInternal" << std::endl;

  char version[SMALL_BUFFER_LENGTH] = "";
  char chromosome[SMALL_BUFFER_LENGTH] = "";
  char line[MED_BUFFER_LENGTH] = "";

  char *tempDirectory = tempnam(NULL, "sw2m.");
  if(mkdir(tempDirectory, S_IRWXU)) {
    free(tempDirectory);
    PROBLEM "Could not create temporary directory %s\nsw2mloci.cpp key 122",
      tempDirectory RECOVER(NULL_ENTRY);
  }
   
  /**
   * Get the chromosome number and SimWalk version by opening the first of
   * the IBD files, reading for the chromosome number and SimWalk version,
   * and closing the file.
   */
  std::ifstream swIbdFile(swIbdFileNames[0]);
  if(swIbdFile.fail()) {
    PROBLEM "\n%s could not be opened for reading.\nsw2mloci.cpp key 80\n",
      swIbdFileNames[0] RECOVER(NULL_ENTRY);
  }

  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(strstr(line, "This run has the integer label:")) {
      multic_strncpy(chromosome, getLastWord(line), SMALL_BUFFER_LENGTH-1);
      //      cout << "chromosome = " << chromosome << std::endl;
    } else if(strstr(line, "Results of the Identity By Descent Analysis from SimWalk2")) {
      multic_strncpy(version, getLastWord(line), SMALL_BUFFER_LENGTH-1);
      //      cout << "version = " << version << std::endl;
    }

    // Stop reading when necessary.
    if(strcmp(chromosome, "") != 0 && strcmp(version, "") != 0) {
      break;
    }
  }

  swIbdFile.close();

  /**
   * Create a mapping between the cM's found in the IBD files to the name of
   * the mibd file that will hold data at that cM.  Then use that map to
   * create a similar map bewtween the cM's and the open file handle.
   */
  CmFileNameMap cmFileNameMap = getCmFileNameMap(swIbdFileNames[0],
						 *swMapFileName,
						 chromosome);
  CmFileMap cmFileMap = getCmFileMap(cmFileNameMap, tempDirectory);

  // Process all of the files specified in swIbdFileNames.
  for(int i = 0; i < *ibdFileCount; i++) {
    const char *swIbdFileName = swIbdFileNames[i];
    Rcout << swIbdFileName << std::endl; // uncommented by PV to debug
    Rcout << ". ";
    Rcout.flush();
    if( (i+1) % 35 == 0) {
      Rcout << std::endl;
    }

    if(strcmp(version, "2.80") == 0) {
      ibd2locis280(swIbdFileName, cmFileMap);
    } else if(strcmp(version, "2.91") == 0) {
      ibd2locis291(swIbdFileName, cmFileMap);
    } else {      
      for_each(cmFileNameMap.begin(), cmFileNameMap.end(),
	       deleteCmFileNameEntry);
      for_each(cmFileMap.begin(), cmFileMap.end(), deleteCmFileEntry);
      PROBLEM "SimWalk version %s not implemented yet.\nsw2mloci.cpp key 165",
	version RECOVER(NULL_ENTRY);
    }
  }
  Rcout << std::endl;

  // Close std::ofstreams and delete the std::ofstream objects.
  for_each(cmFileMap.begin(), cmFileMap.end(), deleteCmFileEntry);

  // Concatenate all the files into one mloci.out
  lociFiles2mloci(cmFileNameMap, tempDirectory);
  rmdir(tempDirectory);
  free(tempDirectory);

  // Free the space allocated for the mibd file names
  for_each(cmFileNameMap.begin(), cmFileNameMap.end(),
	   deleteCmFileNameEntry);
}

/**
 * ibd2locis280 reads a given SimWalk (version 2.80) IBD file and creates
 * many SOLAR-formatted mibd.* files.  It verifies that each IBD file was
 * made by SimWalk 2.80.
 */
void ibd2locis280(const char *swIbdFileName, CmFileMap &cmFileMap) {
  char pedigree[SMALL_BUFFER_LENGTH] = "";
  char version[SMALL_BUFFER_LENGTH] = "";
  char line[MED_BUFFER_LENGTH] = "";

  std::ifstream swIbdFile(swIbdFileName);
  if(swIbdFile.fail()) {
    PROBLEM "\n%s could not be opened for reading.\nsw2mloci.cpp key 75",
      swIbdFileName RECOVER(NULL_ENTRY);
  }

  // Get the pedigree number and test if the version number is correct.
  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(strstr(line, "The results for the pedigree named:")) {
      multic_strncpy(pedigree, getLastWord(line), SMALL_BUFFER_LENGTH-1);
    } else if(strstr(line, "Results of the Identity By Descent Analysis from SimWalk2")) {
      multic_strncpy(version, getLastWord(line), SMALL_BUFFER_LENGTH-1);
      if(strcmp(version, "2.80") != 0) {
	swIbdFile.close();
	PROBLEM "\n%s is version %s, it was expected to be %s\n",
	  swIbdFileName, version, "2.80" RECOVER(NULL_ENTRY);
      }
    }

    // Stop reading when necessary.
    if(strcmp(pedigree, "") != 0 && strcmp(version, "") != 0) {
      break;
    }    
  }

  // These variables are defined outside the while loop for efficiency only.
  char firstId[MED_BUFFER_LENGTH] = "";
  char secondId[MED_BUFFER_LENGTH] = "";
  char marker[MED_BUFFER_LENGTH] = "";
  char cM[MED_BUFFER_LENGTH] = "";

  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(matchesNewPair280(line)) {
      multic_strncpy(firstId, strtok(line, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
      multic_strncpy(secondId, strtok(NULL, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
      multic_strncpy(marker, strtok(NULL, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
      multic_strncpy(cM, strtok(NULL, TOKEN_DELIMITER), MED_BUFFER_LENGTH-1);

      for(int i = 0; i < 3; i++) {
	strtok(NULL, TOKEN_DELIMITER);
      }
      char *kinshipCoeff = strtok(NULL, TOKEN_DELIMITER);

      std::ofstream *mibd = getMibd(cmFileMap, atof(cM));
      *mibd << pedigree << "-" << firstId << "\t"
	    << pedigree << "-" << secondId << "\t" 
	    << kinshipCoeff << "\t" << -1 << std::endl;

      // PV Here's where the -1's get written

    } else if(matchesNewMarker280(line)) {
      multic_strncpy(marker, strtok(line, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
      multic_strncpy(cM, strtok(NULL, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);

      for(int i = 0; i < 3; i++) {
	strtok(NULL, TOKEN_DELIMITER);
      }
      char *kinshipCoeff = strtok(NULL, TOKEN_DELIMITER);

      std::ofstream *mibd = getMibd(cmFileMap, atof(cM));
      *mibd << pedigree << "-" << firstId << "\t"
	    << pedigree << "-" << secondId << "\t" 
	    << kinshipCoeff << "\t" << -1 << std::endl;

      // PV Here's where the -1's get written

    } else if(matchesNewCM280(line)) {
      multic_strncpy(cM, strtok(line, TOKEN_DELIMITER), MED_BUFFER_LENGTH-1);

      for(int i = 0; i < 3; i++) {
	strtok(NULL, TOKEN_DELIMITER);
      }
      char *kinshipCoeff = strtok(NULL, TOKEN_DELIMITER);

      std::ofstream *mibd = getMibd(cmFileMap, atof(cM));
      *mibd << pedigree << "-" << firstId << "\t"
	    << pedigree << "-" << secondId << "\t" 
	    << kinshipCoeff << "\t" << -1 << std::endl;

      // PV Here's where the -1's get written

    }
  }

  swIbdFile.close();
}

/**
 * ibd2locis291 reads a given SimWalk (version 2.91) IBD file and creates
 * many SOLAR-formatted mibd.* files.  It verifies that each IBD file was
 * made by SimWalk 2.91.
 */
void ibd2locis291(const char *swIbdFileName, CmFileMap &cmFileMap) {
  char pedigree[SMALL_BUFFER_LENGTH] = "";
  char version[SMALL_BUFFER_LENGTH] = "";
  char line[MED_BUFFER_LENGTH] = "";

  std::ifstream swIbdFile(swIbdFileName);
  if(swIbdFile.fail()) {
    PROBLEM "\n%s could not be opened for reading.\nsw2mloci.cpp key 75",
      swIbdFileName RECOVER(NULL_ENTRY);
  }

  // Get the pedigree number and verify that the version number is correct.
  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(strstr(line, "The results for the pedigree named:")) {
      multic_strncpy(pedigree, getLastWord(line), SMALL_BUFFER_LENGTH-1);
    } else if(strstr(line, "Results of the Identity By Descent Analysis from SimWalk2")) {
      multic_strncpy(version, getLastWord(line), SMALL_BUFFER_LENGTH-1);
      if(strcmp(version, "2.91") != 0) {
	swIbdFile.close();
	PROBLEM "\n%s is version %s, it was expected to be %s\nsw2mloci.cpp key 292",
	  swIbdFileName, version, "2.91" RECOVER(NULL_ENTRY);
      }
    }

    // Stop reading when necessary.
    if(strcmp(pedigree, "") != 0 && strcmp(version, "") != 0) {
      break;
    }    
  }

  // These variables are defined outside the while loop for efficiency only.
  char firstId[MED_BUFFER_LENGTH] = "";
  char secondId[MED_BUFFER_LENGTH] = "";
  char marker[MED_BUFFER_LENGTH] = "";
  char cM[MED_BUFFER_LENGTH] = "";

  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(matchesNewPair291(line)) {
      multic_strncpy(firstId, strtok(line, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
      multic_strncpy(secondId, strtok(NULL, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
    } else if(matchesNewMarker291(line)) {
      multic_strncpy(marker, strtok(line, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);
      multic_strncpy(cM, strtok(NULL, TOKEN_DELIMITER),
		     MED_BUFFER_LENGTH-1);

      for(int i = 0; i < 3; i++) {
	strtok(NULL, TOKEN_DELIMITER);
      }
      char *kinshipCoeff = strtok(NULL, TOKEN_DELIMITER);

      std::ofstream *mibd = getMibd(cmFileMap, atof(cM));
      *mibd << pedigree << "-" << firstId << "\t"
	    << pedigree << "-" << secondId << "\t" 
	    << kinshipCoeff << "\t" << -1 << std::endl;

    } else if(matchesNewCM291(line)) {
      strtok(line, TOKEN_DELIMITER);
      multic_strncpy(cM, strtok(NULL, TOKEN_DELIMITER), MED_BUFFER_LENGTH-1);

      for(int i = 0; i < 3; i++) {
	strtok(NULL, TOKEN_DELIMITER);
      }
      char *kinshipCoeff = strtok(NULL, TOKEN_DELIMITER);

      std::ofstream *mibd = getMibd(cmFileMap, atof(cM));
      *mibd << pedigree << "-" << firstId << "\t"
	    << pedigree << "-" << secondId << "\t" 
	// The -1 could be replaced by the family structure kinship coeff
	// from the 2.91 IBD files
	    << kinshipCoeff << "\t" << -1 << std::endl;
    }
  }

  swIbdFile.close();
}

/**
 * getLastWord returns the pointer to the last space-separated string of
 * characters from (typically) a line of text.  This is done using the strtok
 * function.  Thus the argument line IS modified during this function.  It
 * does NOT create new memory for such a value.
 */
char *getLastWord(char *line) {
  char *lastWord = NULL;
  char *tmpLastWord = NULL;

  //  localLine = (char *) malloc((strlen(line) + 1) * sizeof(char));
  //  strcpy(localLine, line);

  lastWord = strtok(line, TOKEN_DELIMITER);
  while( (tmpLastWord = strtok(NULL, TOKEN_DELIMITER)) ) {
    lastWord = tmpLastWord;
  }
  
  return lastWord;
}

/**
 * matchesNewPair280 returns true if a line of text matches:
 * id-1 id-2 marker-name cM unused unused unused kinship-coefficient
 * 7500190   7500183   TA178F11   0.000   1.000   0.000   0.000   0.000
 * false otherwise.
 */
bool matchesNewPair280(const char *line) {
  //  cout << "in matchesNewPair280" << std::endl;
  char *localLine = NULL;
  char *token = NULL;
  int tokenCount = 0;

  if(isEmpty(line)) {
    return false;
  }

  localLine = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(localLine, line);

  // Match the first three alphanumeric words
  for(int i = 0; i < 3; i++) {
    if(i == 0) {
      token = strtok(localLine, TOKEN_DELIMITER);
    } else {
      token = strtok(NULL, TOKEN_DELIMITER);
    }
    tokenCount++;
    if(!isWordAlnum(token)) {
      free(localLine);
      return false;
    }
  }

  // Match the next five decimal words
  for(int i = 0; i < 5; i++) {
    token = strtok(NULL, TOKEN_DELIMITER);
    tokenCount++;
    if(!isWordDecimal(token)) {
      free(localLine);
      return false;
    }
  }

  // Make sure we count all the tokens
  while(strtok(NULL, TOKEN_DELIMITER)) {
    tokenCount++;
  }
  if(tokenCount != 8) {
    free(localLine);
    return false;
  }

  free(localLine);
  return true;
}

/**
 * matchesNewMarker280 returns true if a line of text matches:
 * marker-name cM unused unused unused kinship-coefficient
 * D18S481     4.268      1.000      0.000      0.000      0.000
 * false otherwise.
 */
bool matchesNewMarker280(const char *line) {
  //  cout << "in matchesNewMarker280" << std::endl;
  char *localLine = NULL;
  char *token = NULL;
  int tokenCount = 0;

  if(isEmpty(line)) {
    return false;
  }

  localLine = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(localLine, line);

  // Match the first alphanumeric word
  token = strtok(localLine, TOKEN_DELIMITER);
  tokenCount++;
  if(!isWordAlnum(token)) {
    free(localLine);
    return false;
  }

  // Match the next five decimal words
  for(int i = 0; i < 5; i++) {
    token = strtok(NULL, TOKEN_DELIMITER);
    tokenCount++;
    if(!isWordDecimal(token)) {
      free(localLine);
      return false;
    }
  }

  // Make sure we count all the tokens
  while(strtok(NULL, TOKEN_DELIMITER)) {
    tokenCount++;
  }
  if(tokenCount != 6) {
    free(localLine);
    return false;
  }

  free(localLine);
  return true;
}

/**
 * matchesNewCM280 returns true if a line of text matches:
 * cM unused unused unused kinship-coefficient
 * 1.000      1.000      0.000      0.000      0.000
 * false otherwise.
 */
bool matchesNewCM280(const char *line) {
  //  cout << "in matchesNewCM280" << std::endl;
  char *localLine = NULL;
  char *token = NULL;
  int tokenCount = 0;

  if(isEmpty(line)) {
    return false;
  }

  localLine = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(localLine, line);

  // Match the next five decimal words
  for(int i = 0; i < 5; i++) {
    token = strtok(localLine, TOKEN_DELIMITER);
    tokenCount++;
    if(!isWordDecimal(token)) {
      free(localLine);
      return false;
    }
  }

  // Make sure we count all the tokens
  while(strtok(NULL, TOKEN_DELIMITER)) {
    tokenCount++;
  }
  if(tokenCount != 5) {
    free(localLine);
    return false;
  }

  free(localLine);
  return true;
}

/**
 * isWordAlnum returns true if every character in a string passes the alnum
 * test (alpha-numeric), false otherwise.
 */
bool isWordAlnum(const char *word) {
  if(isEmpty(word)) {
    return false;
  }

  for(unsigned int i = 0; i < strlen(word); i++) {
    if(!isalnum(word[i])) {
      return false;
    }
  }
  return true;
}

/**
 * isWordDecimal returns true if every character in a string passes the
 * isdigit test.  word is allowed one decimal point (.).  isWordDecimal
 * returns false if it doesn't fit the above allowances or is the empty
 * string.
 */
bool isWordDecimal(const char *word) {
  if(isEmpty(word)) {
    return false;
  }

  bool foundDecimal = false;

  for(unsigned int i = 0; i < strlen(word); i++) {
    if(!isdigit(word[i])) {
      if(word[i] == '.') {
	if(foundDecimal) {
	  return false;
	}
	foundDecimal = true;
      } else {
	return false;
      }
    }
  }
  return true;
}

/**
 * isEmpty returns true if word is made up of only spaces, tabs, and
 * newlines.
 */
bool isEmpty(const char *word) {
  for(unsigned int i = 0; i < strlen(word); i++) {
    if(!isspace(word[i])) {
      return false;
    }
  }
  return true;
}

/**
 * getMibd is a simple function to perform the find on the CmFileMap and its
 * associated error checks.
 */
std::ofstream *getMibd(CmFileMap &cmFileMap, double cM) {
  CmFileMap::iterator iter = cmFileMap.find(cM);
  if(iter == cmFileMap.end()) {
    PROBLEM "Cannot find key: %f in cM-file map.\nsw2mloci.cpp key 499", cM
      RECOVER(NULL_ENTRY);
  }

  std::ofstream *mibd = iter->second;

  return mibd;
}

/**
 * getSwMibdMap returns a map between the centimorgan (i.e. cM or Kosambi)
 * values in the IBD files to the mibd file that will hold the data for that
 * location.
 *
 * It now detetermines which version of SimWalk make the IBD file and calls
 * a specialized function to parse it and generate the mapping.
 * 
 * Here is the original design:
 * If no map,
 *   - Read first section of ibdFile and simply use 0-x (int) for each
 *     calculated cM.
 * Don't worry about the file existing or not, this is determined in S
 * If .map is good and exists
 *   - Read .map and make a map of <name, kosambe>
 *   - Read first section of ibdFile, counting cM's between markers, keeping
 *     track of the cM's in a vector<double>
 *     - When a new marker appears, add to the map <ibd cM, calculated cM>
 *       and clear vector<double>
 * Return map <ibd cM - double, calculated cM - string>
 */
SwMibdMap getSwMibdMap(const char *mapFileName, const char *ibdFileName) {

  std::ifstream swIbdFile(ibdFileName);
  if(swIbdFile.fail()) {
    PROBLEM "\nCannot open %s for reading.\nsw2mloci.cpp key 614\n",
      ibdFileName RECOVER(NULL_ENTRY);    
  }
  char line[MED_BUFFER_LENGTH];

  SwMibdMap swMibdMap;

  char version[SMALL_BUFFER_LENGTH] = "";
  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(strstr(line, "Results of the Identity By Descent Analysis from SimWalk2")) {
      multic_strncpy(version, getLastWord(line), SMALL_BUFFER_LENGTH-1);
      break;
    }
  }

  if(strcmp(version, "2.80") == 0) {
    swMibdMap = getSwMibdMap280(swIbdFile, mapFileName);
  } else if(strcmp(version, "2.91") == 0) {
    swMibdMap = getSwMibdMap291(swIbdFile, mapFileName);
  } else {
    PROBLEM "SimWalk version %s not implemented yet.\nsw2mloci.cpp key 641",
      version RECOVER(NULL_ENTRY);
  }

  swIbdFile.close();

  /*
    cout << "swMibdMap = " << std::endl;
    for_each(swMibdMap.begin(), swMibdMap.end(), printSwMibdEntry);
  */

  return swMibdMap;
}

/**
 * getSwMibdMap280 specifically reads a SimWalk 2.80 version IBD file and 
 * returns an SwMibdMap
 */
SwMibdMap getSwMibdMap280(std::ifstream &swIbdFile, const char *mapFileName) {
  SwMibdMap swMibdMap;

  NameKosambiMap nameKosambiMap;
  if(strcmp("", mapFileName) != 0) {
    nameKosambiMap = getNameKosambiMap(mapFileName);
  }
  
  char line[MED_BUFFER_LENGTH];
  bool foundFirstPair = false;

  // cMs is used to keep the cM's from between the previous marker and the
  // next.
  vector<double> cMs;
  // previousName is used to remember the previous marker name (somthing like
  // T12AG87)
  char *previousName = NULL;

  // mibdCounter is used when there is no map file, so the mibd files can
  // be named mibd.*.0, mibd.*.1, etc.
  int mibdCounter = 0;

  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(matchesNewPair280(line)) {
      if(foundFirstPair) {
	break;
      }
      foundFirstPair = true;

      strtok(line, TOKEN_DELIMITER); // firstId
      strtok(NULL, TOKEN_DELIMITER); // secondId
      char *tmpPreviousName = strtok(NULL, TOKEN_DELIMITER);
      previousName = (char *) malloc( (strlen(tmpPreviousName) + 1)
				      * sizeof(char));
      strcpy(previousName, tmpPreviousName);

      char *cM = strtok(NULL, TOKEN_DELIMITER);
      cMs.push_back(atof(cM));
    } else if(matchesNewMarker280(line)) {
      char *currentName = strtok(line, TOKEN_DELIMITER);
      char *cM = strtok(NULL, TOKEN_DELIMITER);
      
      if(strcmp("", mapFileName) == 0) {
	for(unsigned int i = 0; i < cMs.size(); i++) {
	  char *stringCm = int2ascii(mibdCounter++);
	  swMibdMap.insert(SwMibdMapEntry(cMs[i], stringCm));
	}
      } else {
	addCmsToMap(swMibdMap, nameKosambiMap, cMs, previousName,
		    currentName);

	free(previousName);      
	previousName = (char *) malloc( (strlen(currentName) + 1)
					* sizeof(char));
	strcpy(previousName, currentName);
      }
     
      cMs.clear();
      cMs.push_back(atof(cM));
    } else if(matchesNewCM280(line)) {
      char *cM = strtok(line, TOKEN_DELIMITER);
      
      cMs.push_back(atof(cM));
    }
  }

  if(strcmp("", mapFileName) == 0) {
    char *stringCm = int2ascii(mibdCounter++);
    swMibdMap.insert(SwMibdMapEntry(cMs[0], stringCm));
  } else {
    addLastCmToMap(swMibdMap, nameKosambiMap, cMs[0], previousName);

    for_each(nameKosambiMap.begin(), nameKosambiMap.end(), 
	     deleteNameKosambiEntry);
  }
  free(previousName);

  /*
    cout << "swMibdMap = " << std::endl;
    for_each(swMibdMap.begin(), swMibdMap.end(), printSwMibdEntry);
  */

  return swMibdMap;
}


/**
 * printSwMibdEntry is used to disply the values in a SwMibdEntry, or all
 * SwMibdEntry's via a for_each command;
 */
void printSwMibdEntry(SwMibdMapEntry &ent) {
  Rcout << "orig cM = " << ent.first << ", "
    "calc cM = " << ent.second << std::endl;  
}

/**
 * deleteSwMibdEntry is used to free the char * memory used for the * in the
 * mibd file name mibd.__.*
 */
void deleteSwMibdEntry(SwMibdMapEntry &ent) {
  free(ent.second);
}

/**
 * getNameKosambiMap creates a map object between the marker name (TA178F11),
 * string, and the centimorgan value it represents, double.  It assumes that
 * the file is a valid .map file.
 */
NameKosambiMap getNameKosambiMap(const char *mapFileName) {
  NameKosambiMap nameKosambiMap;

  std::ifstream mapFile(mapFileName);
  if(mapFile.fail()) {
    PROBLEM "\n%s could not be opened for reading.\nsw2mloci.cpp key 378",
      mapFileName RECOVER(NULL_ENTRY);    
  }

  // Ignore the first line: CHROMOSOME KOSAMBI NAME
  char line[MED_BUFFER_LENGTH];
  mapFile.getline(line, MED_BUFFER_LENGTH);
  
  while(mapFile.getline(line, MED_BUFFER_LENGTH)) {
    if(strcmp(line, "") == 0) {
      break;
    }

    strtok(line, " \t"); // chromosome
    char *kosambi = strtok(NULL, " \t");
    char *name = strtok(NULL, " \t");
    if(strlen(name) > 8) {
      name = &name[strlen(name) - 8];
    }

    char *nameCopy = (char *) malloc((strlen(name) + 1) * sizeof(char));
    strcpy(nameCopy, name);  
  
    nameKosambiMap.insert(NameKosambiMapEntry(nameCopy, atof(kosambi)));
  }

  /*
    cout << "nameKosambiMap = " << std::endl;
    for_each(nameKosambiMap.begin(), nameKosambiMap.end(),
    printNameKosambiEntry);
  */

  return nameKosambiMap;
}

/**
 * printNameKosambiEntry is used to disply the values in a NameKosambiEntry,
 * or all NameKosambiEntry's via a for_each command;
 */
void printNameKosambiEntry(NameKosambiMapEntry &ent) {
  Rcout << "name = " << ent.first << ", "
    "kosambi = " << ent.second << std::endl;  
}

void deleteNameKosambiEntry(NameKosambiMapEntry &ent) {
  //free(ent.first);
}

/**
 * trimTrailingDecimalZeros trims all trailing 0's behind a decimal
 * point (.).  If all digits are trim'ed, it keeps the last 0.
 * Ex: 1.01230 -> 1.0123
 * Ex: 1.0000 -> 1.0
 * Ex: 1 -> 1
 * Ex: 1. -> 1.0
 */
void trimTrailingDecimalZeros(char *decimalNumber) {
  if( !strstr(decimalNumber, ".") ) {
    return;
  }

  while('0' == decimalNumber[strlen(decimalNumber) - 1]) {
    decimalNumber[strlen(decimalNumber) - 1] = '\0';
  }
  if('.' == decimalNumber[strlen(decimalNumber) - 1]) {
    decimalNumber[strlen(decimalNumber)] = '0'; 
    decimalNumber[strlen(decimalNumber)] = '\0'; 
  }
}

/**
 * getCmFileMap creates a map object between the centimorgan value, double,
 * and the open file handle to the desired file, std::ofstream *.  Whether the
 * file name is of the form mibd.(chromosome).(integer) or
 * mibd.(chromosome).(double), depends on whether or not a map file name has
 * been specified
 */
// This probably can get the chromosome itself.
CmFileMap getCmFileMap(CmFileNameMap &cmFileNameMap, const char *directory) {

  CmFileMap cmFileMap;

  for(CmFileNameMap::iterator iter = cmFileNameMap.begin(); 
      iter != cmFileNameMap.end();
      iter++) {

    char *mibdFileName = iter->second;
    int fullFileNameLength = strlen(mibdFileName) + strlen(directory)  + 2;
    char *fullFileName = (char *) malloc(fullFileNameLength
					 * sizeof(char));
    snprintf(fullFileName, fullFileNameLength, "%s/%s", directory,
	     mibdFileName);

    // This could be a temp file instead of a real mibd.(chrom).(cm) file.
    //    std::ofstream *mibd = new std::ofstream(mibdFileName, ios::app);
    std::ofstream *mibd = new std::ofstream(fullFileName, ios::app);
    if(mibd->fail()) {
      PROBLEM "\n%s could not be opened for writing.\nsw2mloci.cpp key 752",
	fullFileName RECOVER(NULL_ENTRY);
    }
    free(fullFileName);

    *mibd << "# " << mibdFileName << std::endl;

    cmFileMap.insert(CmFileMapEntry(iter->first, mibd));
  }

  /*
    cout << "cmFileMap = " << std::endl;
    for_each(cmFileMap.begin(), cmFileMap.end(), printCmFileEntry);
  */

  return cmFileMap;  
}

/**
 * printCmFileEntry is used to disply the values in a CmFileMapEntry,
 * or all CmFileMapEntry's via a for_each command.
 */
void printCmFileEntry(CmFileMapEntry &ent) {
  Rcout << "cM = " << ent.first << ", mibd file address = " << ent.second
       << std::endl;
}

/**
 * deleteCmFileEntry is used to close the mibd file and free the object's
 * memory.
 */
void deleteCmFileEntry(CmFileMapEntry &ent) {
  ent.second->close();
  delete ent.second;
}

/**
 * Concatenate all the files into one mloci.out
 */
void lociFiles2mloci(CmFileNameMap &cmFileNameMap, const char *directory) {
  std::ofstream mloci("mloci.out");
  if(mloci.fail()) {
    PROBLEM "\nError opening mloci.out for writing.\nsw2mloci.cpp key 938"
      RECOVER(NULL_ENTRY);
  }

  for(CmFileNameMap::iterator iter = cmFileNameMap.begin(); 
      iter != cmFileNameMap.end();
      iter++) {
    int fullFileNameLength = strlen(iter->second) + strlen(directory) + 2;
    char *fullFileName = (char *) malloc(fullFileNameLength
					 * sizeof(char));
    snprintf(fullFileName, fullFileNameLength, "%s/%s", directory,
	     iter->second);

    std::ifstream mibd(fullFileName);
    if(mibd.fail()) {
      PROBLEM "\nError opening %s for reading.\nsw2mloci.cpp key 953",
	fullFileName RECOVER(NULL_ENTRY);
    }

    char buffer[MED_BUFFER_LENGTH];
    while(mibd.good()) {
      mibd.getline(buffer, MED_BUFFER_LENGTH - 1);
      mloci << buffer << std::endl;
    }

    mibd.close();
    unlink(fullFileName);
    free(fullFileName);
  }

  mloci.close();
}

// chromosome COULD be fetched from within the function, and not as an
// argument.
CmFileNameMap getCmFileNameMap(const char *swIbdFileName,
			       const char *swMapFileName,
			       const char *chromosome) {
  // First, get a map between the IBD centimorgan values and the appropriate
  // value for the last part of each mibd file name.
  SwMibdMap swMibdMap = getSwMibdMap(swMapFileName, swIbdFileName);

  CmFileNameMap cmFileNameMap;

  for(SwMibdMap::iterator iter = swMibdMap.begin(); 
      iter != swMibdMap.end();
      iter++) {
    char *mibdFileName = (char *) malloc(MED_BUFFER_LENGTH * sizeof(char));
    snprintf(mibdFileName, MED_BUFFER_LENGTH-1, "mibd.%s.%s", chromosome,
	     iter->second);

    cmFileNameMap.insert(CmFileNameMapEntry(iter->first, mibdFileName));
  }

  /*
    cout << "cmFileNameMap = " << std::endl;
    for_each(cmFileNameMap.begin(), cmFileNameMap.end(),
    printCmFileNameEntry);
  */

  return cmFileNameMap;  
}

/**
 * printCmFileNameEntry is used to disply the values in a CmFilenameMapEntry,
 * or all CmFileNameMapEntry's via a for_each command.
 */
void printCmFileNameEntry(CmFileNameMapEntry &ent) {
  Rcout << "cM = " << ent.first << ", mibd file name = " << ent.second
       << std::endl;
}

/**
 * deleteCmFileNameEntry is used to free the char * memory used for the
 * mibd file name.
 */
void deleteCmFileNameEntry(CmFileNameMapEntry &ent) {
  free(ent.second);
}

/**
 * getSwMibdMap291 specifically reads a SimWalk 2.91 version IBD file and 
 * returns an SwMibdMap
 */
SwMibdMap getSwMibdMap291(std::ifstream &swIbdFile, const char *mapFileName) {
  SwMibdMap swMibdMap;

  NameKosambiMap nameKosambiMap;
  if(strcmp("", mapFileName) != 0) {
    nameKosambiMap = getNameKosambiMap(mapFileName);
  }
  
  char line[MED_BUFFER_LENGTH];
  bool foundFirstPair = false;

  // cMs is used to keep the cM's from between the previous marker and the
  // next.
  vector<double> cMs;
  // previousName is used to remember the previous marker found in the IBD
  char *previousName = NULL;

  // mibdCounter is used when there is no map file, so the mibd files can
  // be named mibd.*.0, mibd.*.1, etc.
  int mibdCounter = 0;

  while(swIbdFile.getline(line, MED_BUFFER_LENGTH)) {
    if(matchesNewPair291(line)) {
      if(foundFirstPair) {
	break;
      }
      foundFirstPair = true;
    } else if(matchesNewMarker291(line)) {
      if( !previousName ) {
	char *tmpPreviousName = strtok(line, TOKEN_DELIMITER);
	previousName = (char *) malloc( (strlen(tmpPreviousName) + 1)
					* sizeof(char));
	strcpy(previousName, tmpPreviousName);
	
	char *cM = strtok(NULL, TOKEN_DELIMITER);
	cMs.push_back(atof(cM));
	
	continue;
      }
      char *currentName = strtok(line, TOKEN_DELIMITER);
      char *cM = strtok(NULL, TOKEN_DELIMITER);
      
      if(strcmp("", mapFileName) == 0) {
	for(unsigned int i = 0; i < cMs.size(); i++) {
	  char *stringCm = int2ascii(mibdCounter++);
	  swMibdMap.insert(SwMibdMapEntry(cMs[i], stringCm));
	}
      } else {
	addCmsToMap(swMibdMap, nameKosambiMap, cMs, previousName,
		    currentName);
	
	free(previousName);
	previousName = (char *) malloc( (strlen(currentName) + 1)
					* sizeof(char));
	strcpy(previousName, currentName);
      }
      
      cMs.clear();
      cMs.push_back(atof(cM));
    } else if(matchesNewCM291(line)) {
      strtok(line, TOKEN_DELIMITER); // dashes
      char *cM = strtok(NULL, TOKEN_DELIMITER);
      
      cMs.push_back(atof(cM));
   }
  }

  if(strcmp("", mapFileName) == 0) {
    char *stringCm = int2ascii(mibdCounter++);
    swMibdMap.insert(SwMibdMapEntry(cMs[0], stringCm));
  } else {
    addLastCmToMap(swMibdMap, nameKosambiMap, cMs[0], previousName);
  
    for_each(nameKosambiMap.begin(), nameKosambiMap.end(), 
	     deleteNameKosambiEntry);
  }
  free(previousName);

  return swMibdMap;
}

/**
 * matchesNewPair291 returns true if a line of text matches:
 * id-1 id-2 kinship-coefficient
 * 7500190   7500183   0.000
 * false otherwise.
 */
bool matchesNewPair291(const char *line) {
  //  cout << "in matchesNewPair291" << std::endl;
  char *localLine = NULL;
  char *token = NULL;
  int tokenCount = 0;

  if(isEmpty(line)) {
    return false;
  }

  localLine = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(localLine, line);

  // Match the first two alphanumeric words
  for(int i = 0; i < 2; i++) {
    if(i == 0) {
      token = strtok(localLine, TOKEN_DELIMITER);
    } else {
      token = strtok(NULL, TOKEN_DELIMITER);
    }
    tokenCount++;
    if(!isWordAlnum(token)) {
      free(localLine);
      return false;
    }
  }

  // Match the next decimal word
  token = strtok(NULL, TOKEN_DELIMITER);
  tokenCount++;
  if(!isWordDecimal(token)) {
    free(localLine);
    return false;
  }

  // Make sure we count all the tokens
  while(strtok(NULL, TOKEN_DELIMITER)) {
    tokenCount++;
  }
  if(tokenCount != 3) {
    free(localLine);
    return false;
  }

  free(localLine);
  return true;
}

/**
 * matchesNewMarker291 returns true if a line of text matches:
 * marker-name cM unused unused unused kinship-coefficient
 * D18S481     4.268      1.000      0.000      0.000      0.000
 * false otherwise.
 *
 * As of 2006-01-23, there are no differences between version 2.80 and 2.91
 * SimWalk IBDs as how a new marker is formatted.  If this changes, this
 * function will need to not call the 2.80 version.
 */
bool matchesNewMarker291(const char *line) {
  //  cout << "in matchesNewMarker291" << std::endl;
  return matchesNewMarker280(line);
}

/**
 * matchesNewCM291 returns true if a line of text matches:
 * cM unused unused unused kinship-coefficient
 * --     1.000      1.000      0.000      0.000      0.000
 * false otherwise.
 */
bool matchesNewCM291(const char *line) {
  //  cout << "in matchesNewCM291" << std::endl;
  char *localLine = NULL;
  char *token = NULL;
  int tokenCount = 0;

  if(isEmpty(line)) {
    return false;
  }

  localLine = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(localLine, line);

  // In version 2.91, each new cM begins with a --
  if(strcmp(strtok(localLine, TOKEN_DELIMITER), "--") != 0) {    
    return false;
  }
  tokenCount++;

  // Match the next five decimal words
  for(int i = 0; i < 5; i++) {
    token = strtok(NULL, TOKEN_DELIMITER);
    tokenCount++;
    if(!isWordDecimal(token)) {
      free(localLine);
      return false;
    }
  }

  // Make sure we count all the tokens
  while(strtok(NULL, TOKEN_DELIMITER)) {
    tokenCount++;
  }
  if(tokenCount != 6) {
    free(localLine);
    return false;
  }

  free(localLine);
  return true;
}

/**
 * addCmsToMap is used to aggragate the steps needed to add the
 * SimWalk cM-SOLAR cM pairs to the SwMibdMap.  This is done in all versions
 * of getSwMibdMap*, so it was better to pull it out.
 */
void addCmsToMap(SwMibdMap &swMibdMap, NameKosambiMap &nameKosambiMap,
		 vector<double> &cMs,		 
		 const char *previousName, const char *currentName) {
  NameKosambiMap::iterator iter = nameKosambiMap.find(previousName);
  if(iter == nameKosambiMap.end()) {
    PROBLEM "\nCannot find key: %s in name-kosambi map.\nsw2mloci.cpp key 449",
      previousName RECOVER(NULL_ENTRY);
  }
  double previousCm = iter->second;
  
  iter = nameKosambiMap.find(currentName);
  if(iter == nameKosambiMap.end()) {
    PROBLEM "\nCannot find key: %s in name-kosambi map.\n sw2mloci.cpp key 546",
      currentName RECOVER(NULL_ENTRY);
  }
  double currentCm = iter->second;
  
  if(currentCm <= previousCm) {
    PROBLEM "\nThe key: %s (%f) does not have a larger centimorgan value\nthan the key: %s (%f)\nThe map file is not formatted properly.\nsw2mloci.cpp key 552.",
      currentName, currentCm, previousName, previousCm RECOVER(NULL_ENTRY);
  }
  
  double delta = (currentCm - previousCm) / cMs.size();
  
  for(unsigned int i = 0; i < cMs.size(); i++) {
    char *stringCm = double2ascii(previousCm + i * delta, 3);
    swMibdMap.insert(SwMibdMapEntry(cMs[i], stringCm));
  }
}

/**
 * addLastCmToMap is used to aggragate the steps needed to add the last
 * marker data to the SwMibdMap.  This is done in all versions of
 * getSwMibdMap*, so it was better to pull it out.
 */
void addLastCmToMap(SwMibdMap &swMibdMap, NameKosambiMap &nameKosambiMap,
		    double cM, const char *markerName) {
  NameKosambiMap::iterator iter = nameKosambiMap.find(markerName);
  if(iter == nameKosambiMap.end()) {
    PROBLEM "%s is not a key in the name-kosambi map.\nsw2mloci.cpp key 1299",
      markerName RECOVER(NULL_ENTRY);
  }
  double markerCm = iter->second;
  
  char *stringCm = double2ascii(markerCm, 3);
  
  swMibdMap.insert(SwMibdMapEntry(cM, stringCm));
}

/**
 * int2ascii aggragates the computations to transform an int into a string.
 */
char *int2ascii(int i) {
  char tmpString[SMALL_BUFFER_LENGTH];
  snprintf(tmpString, SMALL_BUFFER_LENGTH-1, "%d", i);

  char *string = (char *) malloc( (strlen(tmpString) + 1) * sizeof(char));
  multic_strncpy(string, tmpString, strlen(tmpString));

  return string;
}

/**
 * double2ascii aggragates the computations to transform an double into a
 * string.  The precision (digits past the decimal) can be specified.
 */
char *double2ascii(double d, int precision) {
  if(precision < 0) {
    PROBLEM "precision (%d) less than 0 is not allowed.\nsw2mloci.cpp key 1219",
      precision RECOVER(NULL_ENTRY);
  }

  char formatString[SMALL_BUFFER_LENGTH];
  snprintf(formatString, SMALL_BUFFER_LENGTH-1, "%%.%df", precision);

  char tmpString[SMALL_BUFFER_LENGTH];
  snprintf(tmpString, SMALL_BUFFER_LENGTH-1, formatString, d);
  trimTrailingDecimalZeros(tmpString);

  char *string = (char *) malloc( (strlen(tmpString) + 1) * sizeof(char));
  multic_strncpy(string, tmpString, strlen(tmpString));
  
  return string;
}

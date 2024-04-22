/***************************************************************************
 *   Copyright (C) 2021 by Marco Beccuti				                   *
 *   marco.beccuti@unito.it						                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef __CGLPK_H__
	#define __CGLPK_H__
	#include "GLPKsolve.hpp"
#endif


namespace FBGLPK{


	/**
	 * @brief Parses flux names from an input stream and updates internal mappings.
	 * 
	 * This method reads a single line from the provided input stream, which is expected to contain flux names separated by a space.
	 * Each flux name is then stored in two internal structures: one mapping the names to their corresponding indices, and another 
	 * maintaining the order of these names as they appear in the input. If the line of flux names cannot be read or if no names are found,
	 * the method throws an exception.
	 *
	 * @param in Reference to an input file stream from which flux names are read.
	 * @param parser A parser object used to split the input line into names based on provided delimiters.
	 * @param delimC A C-style string containing delimiter characters that separate flux names in the input.
	 * 
	 * @throws Exception If the line of flux names cannot be read or if it contains no recognizable names.
	 * 
	 */
	void LPprob::parseFluxNames(std::ifstream& in, general::Parser& parser, const char* delimC) {
		  std::string buffer;
		  if (!std::getline(in, buffer)) {
		      throw Exception("FLUX BALANCE: Failed to read flux names line.");
		  }
		  parser.update(delimC, buffer);
		  if (parser.size() == 0) {
		      throw Exception("FLUX BALANCE: No flux names found.");
		  }
		  for (unsigned int i = 0; i < parser.size(); ++i) {
		      ReactionsNamesId[parser.get(i)] = i + 1;
		      ReactionsNamesOrd.push_back(parser.get(i));
		  }
	}


	/**
	 * @brief Parses the matrix dimensions and objective type from an input stream.
	 * 
	 * This method reads a line from the provided input stream expected to contain three parts: number of rows, number of columns,
	 * and a string representing the objective type of the optimization problem (GLP_MAX,). 
	 * It then updates the provided sizeRow and sizeCol parameters with these values. 
	 * If variability is not flagged (variability == 0), it also sets the typeOBJ parameter.
	 *
	 * @param in Reference to an input file stream from which the dimensions and type are read.
	 * @param parser A parser object used to split the input line into parts based on the provided delimiters.
	 * @param delimC A C-style string containing delimiter characters that separate the data items in the input.
	 * @param sizeRow Reference to an unsigned int where the number of rows will be stored.
	 * @param sizeCol Reference to an unsigned int where the number of columns will be stored.
	 * @param typeOBJ Reference to an integer where the objective type of the LP problem will be stored.
	 * @param variability An integer flag indicating whether the objective type should be set (0 if it should be set).
	 * 
	 * @throws Exception If the line cannot be read, if it does not contain exactly three parts, or if the parts cannot be
	 * correctly converted and interpreted.
	 *
	 */

	void LPprob::parseModelDimensionsAndType(std::ifstream& in, general::Parser& parser, const char* delimC, unsigned int& sizeRow, unsigned int& sizeCol, int& typeOBJ, int variability) {
		  std::string buffer;
		  if (!std::getline(in, buffer)) {
		      throw Exception("FLUX BALANCE: Failed to read dimensions and type line.");
		  }
		  parser.update(delimC, buffer);
		  if (parser.size() != 3) {
		      throw Exception("FLUX BALANCE: Incorrect format for dimensions and type line.");
		  }
		  sizeRow = std::stoul(parser.get(0));
		  sizeCol = std::stoul(parser.get(1));
		  if(!variability)
		  	typeOBJ = setTypeObj(parser.get(2));
	}


	/**
	 * @brief Parses and sets the objective coefficients for an LP problem from an input stream.
	 *
	 * This method reads a line of objective coefficients from the specified input stream, expecting them to be 
	 * separated by the given delimiters. It then checks if the number of coefficients matches the expected number of columns.
	 * If `setDefaultCoefficients` is true, it directly sets these coefficients in the LP problem. If `variability` is true,
	 * it instead sets a specific coefficient (related to `flux_var`) to 1.0, making the corresponding flux the focus of optimization.
	 *
	 * @param in Reference to the input stream from which coefficients are read.
	 * @param parser Parser object used to split the input string.
	 * @param delimC String of delimiter characters for parsing the input.
	 * @param setDefaultCoefficients Flag indicating if coefficients should be directly set.
	 * @param variability Flag indicating if the setup should adapt based on varying conditions.
	 * @param flux_var Index specifying which flux's coefficient to set to 1.0 under variability.
	 * @param var_obj_eq Optional pointer to a string where the raw input line is stored if variability is active.
	 *
	 * @throws Exception If unable to read the line or if the number of parsed coefficients does not match the expected number.
	 */
	void LPprob::parseObjectiveCoefficients(std::ifstream& in, general::Parser& parser, const char* delimC, 
		                                      int setDefaultCoefficients, int variability, 
		                                      int flux_var, std::string* var_obj_eq) {
		  std::string buffer;
		  if (!std::getline(in, buffer)) {
		      throw Exception("FLUX BALANCE: Failed to read objective coefficients.");
		  }
		  parser.update(delimC, buffer);
		  if (parser.size() != sizeCol) {
		      throw Exception("FLUX BALANCE: Incorrect number of objective coefficients.");
		  }
		  if (setDefaultCoefficients) {
		      for (unsigned int i = 0; i < parser.size(); ++i) {
		          glp_set_obj_coef(lp, i + 1, std::stof(parser.get(i)));
		      }
		  } else if (variability) {
		      if (var_obj_eq != nullptr) {
		          *var_obj_eq = buffer;  // Modifica sicura solo se var_obj_eq non è nullptr
		      }
		      for (unsigned int i = 0; i < sizeCol; ++i) {
		          glp_set_obj_coef(lp, i + 1, (i + 1 == flux_var) ? 1.0 : 0.0);
		      }
		  }
	}

	/**
	 * @brief Sets the bounds for each row (constraint) in the linear programming problem.
	 *
	 * This method reads bounds for each row from an input stream, expecting each line to specify the bound type and limits for one row.
	 * The bounds are defined by three parts: the bound type, lower limit, and upper limit, separated by specified delimiters. 
	 * If the format of any line does not match the expected format (three parts), it throws an exception. 
	 * The bounds types are converted from string identifiers to GLPK-specific constants using the `setTypeBound` method.
	 *
	 * @param in Reference to the input file stream from which row bounds are read.
	 * @param parser Parser object used to split the input line into components.
	 * @param delimC String of delimiter characters for parsing the input.
	 *
	 * @throws Exception If any line does not contain exactly three parts or fails to meet format expectations.
	 *
	 */
	 
	void LPprob::setRowBounds(std::ifstream& in, general::Parser& parser, const char* delimC) {
		  std::string buffer;
		  for (unsigned int i = 0; i < sizeRow && std::getline(in, buffer); ++i) {
		      parser.update(delimC, buffer);
		      if (parser.size() != 3) {
		          throw Exception("FLUX BALANCE: Incorrect row bounds format.");
		      }
		      glp_set_row_bnds(lp, i + 1, setTypeBound(parser.get(0)), std::atof(parser.get(1).c_str()), std::atof(parser.get(2).c_str()));
		  }
	}


	/**
	 * @brief Sets the bounds for each column (variable) in the linear programming problem.
	 *
	 * This method reads bounds for each column from an input stream, with each line expected to specify the bound type, 
	 * lower limit, and upper limit for one column. These parts should be separated by specified delimiters. The method processes 
	 * each line to set the bounds in the LP problem using GLPK functions. If a line does not match the expected format (three parts), 
	 * an exception is thrown. Special handling is implemented for GLP_DB type bounds: if both lower and upper bounds are zero, 
	 * the column is set as fixed (GLP_FX), meaning its value cannot vary.
	 *
	 * @param in Reference to the input file stream from which column bounds are read.
	 * @param parser Parser object used to split the input line into components.
	 * @param delimC String of delimiter characters for parsing the input.
	 *
	 * @throws Exception If any line does not contain exactly three parts or fails to meet format expectations.
	 *
	 */
	void LPprob::setColumnBounds(std::ifstream& in, general::Parser& parser, const char* delimC) {
		  std::string buffer;
		  for (unsigned int i = 0; i < sizeCol && std::getline(in, buffer); ++i) {
		      parser.update(delimC, buffer);
		      if (parser.size() != 3) {
		          throw Exception("FLUX BALANCE: Incorrect column bounds format.");
		      }
		      double lb = std::atof(parser.get(1).c_str());
		      double ub = std::atof(parser.get(2).c_str());
		      if (parser.get(0) == "GLP_DB" && lb == ub) {
		          glp_set_col_bnds(lp, i + 1, GLP_FX, lb, lb);
		      } else {
		          glp_set_col_bnds(lp, i + 1, setTypeBound(parser.get(0)), lb, ub);
		      }
		  }
	}

	/**
	 * @brief Sets the sparse matrix coefficients for the linear programming problem.
	 *
	 * This method reads matrix entries from an input stream, with each line expected to specify the row index, 
	 * column index, and value of a non-zero element in the sparse matrix. These entries are used to construct the 
	 * matrix in the LP problem. If a line does not conform to the expected format (three components), an exception is thrown.
	 * This method populates three arrays: ia, ja, and ar, which represent the row indices, column indices, and values of the
	 * non-zero elements, respectively.
	 *
	 * @param in Reference to the input file stream from which matrix entries are read.
	 * @param parser Parser object used to split the input line into components.
	 * @param delimC String of delimiter characters for parsing the input.
	 *
	 * @throws Exception If any line does not contain exactly three parts or fails to meet format expectations.
	 *
	 */
	void LPprob::setSparseMatrix(std::ifstream& in, general::Parser& parser, const char* delimC) {
		  std::string buffer;
		  for (unsigned int i = 0; i < sizeVet && std::getline(in, buffer); ++i) {
		      parser.update(delimC, buffer);
		      if (parser.size() != 3) {
		          throw Exception("FLUX BALANCE: Incorrect matrix format.");
		      }
		      ia[i + 1] = std::atoi(parser.get(0).c_str());
		      ja[i + 1] = std::atoi(parser.get(1).c_str());
		      ar[i + 1] = std::atof(parser.get(2).c_str());
		  }
	}


	/**
	 * @brief Parses the count of non-zero elements (NZE) from an input stream.
	 * 
	 * This method reads a single line from the provided input stream which should contain the count of non-zero elements 
	 * in the sparse matrix. The count is then converted to an unsigned integer and stored in the sizeVet parameter.
	 *
	 * @param in Reference to an input file stream from which the NZE count is read.
	 * @param parser A parser object used to process the input line and extract the NZE count.
	 * @param delimC A C-style string containing delimiter characters that separate the data items in the input.
	 * @param sizeVet Reference to an unsigned int where the NZE count will be stored.
	 * 
	 * @throws Exception If the NZE count line cannot be read or if the count is not properly formatted.
	 * 
	 * Usage note: This function assumes the NZE count is the only item on the read line, correctly formatted as a numeric string.
	 */

	void LPprob::parseSizeVet(std::ifstream& in, general::Parser& parser, const char* delimC, unsigned int& sizeVet) {
		  std::string buffer;
		  if (!std::getline(in, buffer)) {
		      throw Exception("FLUX BALANCE: Failed to read non-zero elements count.");
		  }
		  parser.update(delimC, buffer);
		  sizeVet = std::stoul(parser.get(0));
	}

	/**
	 * @brief Initializes the linear programming (LP) problem with GLPK.
	 * 
	 * This function sets up a new LP problem instance using the GNU Linear Programming Kit (GLPK). It creates a new problem object
	 * and assigns a name to it.
	 *
	 * @param fileProb The name of the file, which will also be used as the name of the LP problem for identification.
	 * 
	 */
	void LPprob::initializeLP(const char* fileProb) {
		  lp = glp_create_prob();
		  glp_set_prob_name(lp, fileProb);
	}

	/**
	 * @brief Allocates memory for matrix indices and values based on the non-zero elements count.
	 * 
	 * This function dynamically allocates memory for arrays that will store the row indices (ia), column indices (ja),
	 * non-zero element values (ar), and additional values (Value), which are all used in constructing the sparse matrix
	 * for the linear programming problem.
	 * If any allocation fails, all allocated memory is freed, and a std::bad_alloc exception is thrown to signal the error.
	 *
	 * @param sizeVet The number of non-zero elements in the sparse matrix, used to determine the amount of memory to allocate.
	 * 
	 * @throws std::bad_alloc If memory allocation for any of the arrays fails.
	 */
	void LPprob::allocateMemory(unsigned int sizeVet){

		  ia = (int*)malloc(sizeof(int) * (sizeVet + 1));
		  ja = (int*)malloc(sizeof(int) * (sizeVet + 1));
		  ar = (double*)malloc(sizeof(double) * (sizeVet + 1));
		  Value = (double*)malloc(sizeof(double) * (sizeVet + 1));

		  if (!ia || !ja || !ar || !Value) {
		      free(ia); free(ja); free(ar); free(Value);  
		      throw std::bad_alloc(); 
		  }
	}


	/**
	 * @brief Manages the adjustment of objective coefficients based on variability and a specified flux name.
	 *
	 * This method adjusts the linear programming problem's focus based on the presence of variability and a specific flux name.
	 * If variability is enabled, it verifies whether the provided FluxName is a valid flux name present in the internal mapping
	 * of reaction names to indices (ReactionsNamesId). If the flux name is found, it sets the index of this flux as the target
	 * for variability adjustments in the optimization problem. If the flux name is not found, an exception is thrown indicating
	 * the invalid flux name.
	 *
	 * @param in Reference to the input file stream, used if additional information needs to be read (not used in the current implementation).
	 * @param parser A parser object potentially used for further input parsing (not used in the current implementation).
	 * @param delimC A C-style string containing delimiter characters; provided for potential parsing needs (not used in the current implementation).
	 * @param variability An integer indicating if variability adjustments are active (non-zero means active).
	 * @param FluxName A C-string representing the name of the flux to focus on if variability is active.
	 *
	 * @throws Exception If the flux name provided is not found in the internal mappings, indicating an invalid or unspecified flux.
	 *
	 */
	void LPprob::manageVariability(ifstream& in, general::Parser& parser, const char* delimC, int variability, const char* FluxName) {
		  if (variability) {
		      auto it = ReactionsNamesId.find(string(FluxName));
		      if (it == ReactionsNamesId.end())
		          throw Exception(string(FluxName) + " is not a valid flux name");
		      flux_var = it->second;
		  }
	}


	/**
	 * @brief Adds a new row to the sparse matrix to account for variability in the objective function.
	 *
	 * This method integrates additional objective coefficients into the sparse matrix as a new row, 
	 * representing variability adjustments specified in the optimization problem. It reads coefficients
	 * from a given string, checks their count against the expected number of columns, and appends them
	 * as a new row at the end of the existing matrix. This allows for dynamic adjustments to the model
	 * based on specific scenarios or requirements. If the number of coefficients does not match the number 
	 * of columns, an exception is thrown.
	 *
	 * @param var_obj_eq A string containing the coefficients for the variability row, separated by delimiters.
	 * @param parser A parser object used to split the variability coefficients string.
	 * @param delimC A C-style string containing delimiter characters that separate the coefficients.
	 * @param sizeVet Reference to an unsigned integer tracking the total number of non-zero elements in the matrix,
	 *                which will be updated to include the new row's elements.
	 *
	 * @throws Exception If the number of coefficients parsed does not match the expected number of columns.
	 *
	 */
	void LPprob::addVariabilityRow(const std::string& var_obj_eq, general::Parser& parser, const char* delimC, unsigned int& sizeVet) {
		  parser.update(delimC, var_obj_eq);

		  if (parser.size() != sizeCol) {
		      throw Exception("FLUX BALANCE: Incorrect number of coefficients for variability objective.");
		  }

		  int baseIndex = sizeVet + 1;  // Starting index for new entries in ia, ja, ar
		  for (unsigned int j = 0; j < parser.size(); ++j) {
		      int index = baseIndex + j;
		      ia[index] = sizeRow + 1;  // Index for the new variability row
		      ja[index] = j + 1;        // Column index
		      ar[index] = atof(parser.get(j).c_str());
		  }

		  sizeVet += parser.size();  // Update sizeVet to reflect added elements
	}

	/**
	 * @brief Constructs an LPprob object by setting up a linear programming problem.
	 * 
	 * This constructor initializes the linear programming problem by reading from a specified file.
	 * The file should contain all necessary details of the problem including flux names, matrix dimensions,
	 * objective coefficients, and bounds for rows and columns. It sets up the environment, parses and loads
	 * data into the GLPK solver for further operations like optimization. The process involves multiple steps:
	 * opening the file, initializing the problem, parsing the file contents, allocating memory, and setting up 
	 * the GLPK matrix with the provided data.
	 *
	 * @param fileProb Path to the file containing the problem specifications.
	 * 
	 * @exception std::exception Throws if the file cannot be opened, or if any parsing or initialization step fails.
	 * 
	 * Example usage:
	 * @code
	 *   LPprob lpProblem("path/to/problem/file.txt");
	 * @endcode
	 */
	LPprob::LPprob(const char* fileProb){
		  try{
		      ifstream in(fileProb, std::ifstream::in);
		      if (!in) {
		          throw Exception("FLUX BALANCE: error opening input file:" + string(fileProb));
		      }

		      initializeLP(fileProb);

		      general::Parser parser;
		      char delimC[] = "\t, ;\"";
		      int typeOBJ;

		      parseFluxNames(in, parser, delimC);
		      parseModelDimensionsAndType(in, parser, delimC, sizeRow, sizeCol, typeOBJ, 0);
		      parseSizeVet(in, parser, delimC, sizeVet);

		      allocateMemory(sizeVet);

		      glp_set_obj_dir(lp, typeOBJ);
		      glp_add_rows(lp, sizeRow);
		      glp_add_cols(lp, sizeCol);

		      parseObjectiveCoefficients(in, parser, delimC, 1, 0, 0, nullptr);
		      setRowBounds(in, parser, delimC);
		      setColumnBounds(in, parser, delimC);
		      setSparseMatrix(in, parser, delimC);
		      
		      glp_load_matrix(lp, sizeVet, ia, ja, ar);

		  }catch(const std::exception& e){
		      std::cerr << "Exception: " << e.what() << std::endl;
		      exit(EXIT_FAILURE);
		  }
	}
	

	LPprob::LPprob(const char* FileProb, const char* FileInVar, const char* FileOutVar, int typeOBJ, const char* FluxName, const int gamma) {
		  // Creating LP problem
		  updateLP(FileProb, 1, typeOBJ, FluxName);

		  try {
		      // Opening input file
		      in_var.open(FileInVar, std::ifstream::in);
		      if (!in_var)
		          throw Exception("FLUX BALANCE: error opening input file:" + string(FileInVar));
		      if (in_var.eof())
		          throw Exception("FLUX BALANCE: error input file:" + string(FileInVar) + " is empty");

		      this->gamma = gamma;

		      string buffer;
		      getline(in_var, buffer);
		      // Reading flux names
		      // cout << buffer << endl

		      // Opening output file
		      out_var.open(FileOutVar, std::ofstream::out);
		      if (!out_var)
		          throw Exception("FLUX BALANCE: error opening output file:" + string(FileOutVar));
		      out_var.precision(16);
		  } catch (exception& e) {
		      cout << "\nException: " << e.what() << endl;
		      exit(EXIT_FAILURE);
		  } catch (Exception& e) {
		      cout << "\nException: " << e.what() << endl;
		      exit(EXIT_FAILURE);
		  }
	}

 


	/**
	 * @brief Updates or initializes a linear programming problem based on file input, with optional variability adjustments.
	 *
	 * This method configures or updates the settings of an LP problem by reading specifications from a given file.
	 * It involves setting flux names, problem dimensions, objective coefficients, and bounds. The process is sensitive to
	 * the `variability` parameter, which dictates whether additional adjustments are made to cater to dynamic aspects of the
	 * problem setup, such as changing objective coefficients or constraints based on a specific flux name.
	 *
	 * If variability is enabled, an additional row may be added, and specific objective coefficients are adjusted to focus on a designated flux variable.
	 *
	 * @param fileProb The file path from which to read the LP problem settings.
	 * @param variability An integer flag indicating whether to apply variability-specific settings (1 if true).
	 * @param typeOBJ The objective type (maximization or minimization) of the LP problem.
	 * @param FluxName The name of the flux variable to focus on if variability is enabled.
	 *
	 * @throws Exception If the file cannot be opened or parsed correctly, or if any setup step fails due to format issues or
	 *                   logical errors in configuration.
	 *
	 * Usage example:
	 * @code
	 *   LPprob lp;
	 *   lp.updateLP("path/to/problem.txt", 1, GLP_MAX, "targetFlux");
	 * @endcode
	 */
	void LPprob::updateLP(const char* fileProb, int variability, int typeOBJ, const char* FluxName) {
		  try {
		      string var_obj_eq = "";
		      ifstream in(fileProb, std::ifstream::in);
		      if (!in) {
		          throw Exception("FLUX BALANCE: error opening input file:" + string(fileProb));
		      }

		      initializeLP(fileProb);
		      general::Parser parser;
		      char delimC[] = "\t, ;\"";

		      parseFluxNames(in, parser, delimC);
		      parseModelDimensionsAndType(in, parser, delimC, sizeRow, sizeCol, typeOBJ, variability);
		      parseSizeVet(in, parser, delimC, sizeVet);

		      if (variability) {
		          manageVariability(in, parser, delimC, variability, FluxName);
		      }

		      allocateMemory(sizeVet + (variability * sizeCol));
		      glp_set_obj_dir(lp, typeOBJ);
		      glp_add_rows(lp, sizeRow + variability);
		      glp_add_cols(lp, sizeCol);
		      parseObjectiveCoefficients(in, parser, delimC, !variability, variability, flux_var, &var_obj_eq);

		      setRowBounds(in, parser, delimC);
		      setColumnBounds(in, parser, delimC);
		      setSparseMatrix(in, parser, delimC);

		      if (variability) {
		          addVariabilityRow(var_obj_eq, parser, delimC, sizeVet);
		      }

		      glp_load_matrix(lp, sizeVet, ia, ja, ar);

		      filename = string(fileProb);
		      
		  } catch (exception& e) {
		      cout << "Exception: " << e.what() << endl;
		      exit(EXIT_FAILURE);
		  } catch (Exception& e) {
		      cout << "Exception: " << e.what() << endl;
		      exit(EXIT_FAILURE);
		  }
	}



	void LPprob::solveVariability(){
		  string buffer;
		  out_var<<"Time Obj"<<endl;
		  try{
		      while (!in_var.eof()){
		      getline(in_var,buffer); 
		      class general::Parser par(" ",buffer);
		      //cout<<"tot:"<<par.size()<<" sizeCol"<<sizeCol<<endl;
		      if (par.size()!=0){
		       for (unsigned int i=sizeCol+2,j=1; i<par.size();i=i+2,++j){ //+2 is due to time and obj
		          //updating bound
		          update_bound(j,get_bound_type(j),atof(par.get(i).c_str()),atof(par.get(i+1).c_str()));
		          }
		       //updating new equation based on old obj
		       glp_set_row_bnds(lp,sizeRow+1, GLP_LO, atof(par.get(0).c_str())*gamma, atof(par.get(0).c_str())*gamma);
		       //glp_set_col_bnds(lp,sizeRow+1, GLP_LO, atof(par.get(0).c_str())*gamma, atof(par.get(0).c_str())*gamma);
		       solve();
		       out_var<<par.get(0)<<" "<< glp_get_obj_val(lp)<<endl;   
		      }
		      }
		   }
	 catch (exception& e){
		  cout << "\n Exception: " << e.what() << endl;
		  exit(EXIT_FAILURE);
		  }
	 catch (Exception& e){
		  cout << "\n Exception: " << e.what() << endl;
		  exit(EXIT_FAILURE);
		  }       
	}

	int LPprob::setTypeBound(string typeString){
		  int type;
		  if (typeString.compare("GLP_FR")==0)
		      type=GLP_FR;
		  else  if (typeString.compare("GLP_LO")==0)
		      type=GLP_LO;
		  else  if (typeString.compare("GLP_UP")==0)
		      type=GLP_UP;
		  else  if (typeString.compare("GLP_DB")==0)
		      type=GLP_DB;
		  else
		      type=GLP_FX;
		  return type;
	}


}

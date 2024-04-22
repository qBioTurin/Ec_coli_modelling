#ifndef __GLPK__
    #define __GLPK__
    #include <glpk.h>
#endif

#ifndef __IOS_H__
    #define __IOS_H__
    #include <iostream>
#endif

#ifndef __FSTREAM__
    #define __FSTREAM__
    #include <fstream>
#endif

#ifndef __UNORDERED_MAP__
    #define __UNORDERED_MAP__
    #include <unordered_map>
#endif

#ifndef __GEN_H__
    #define __GEN_H__
    #include "general.h"
#endif

namespace FBGLPK {
    using namespace std;

    /**
     * @brief Exception class for handling runtime errors within FBGLPK.
     */
    struct Exception {
        std::string mess; //!< Error message

        Exception() : mess("") {} //!< Default constructor
        Exception(std::string mess) : mess(mess) {} //!< Constructor with error message

        /**
         * @brief Retrieves the error message.
         * @return Error message as a string.
         */
        std::string what() const { return mess; }
    };

    /**
     * @brief Class representing a linear programming problem using GLPK.
     * 
     * This class encapsulates functionalities for setting up and solving linear programming problems,
     * handling variability in solutions, and exporting results. It heavily relies on the GLPK library.
     */
    class LPprob {
    private:
        glp_prob *lp {nullptr}; //!< GLPK problem object
        int *ia {nullptr}; //!< Row indices for the LP matrix
        int *ja {nullptr}; //!< Column indices for the LP matrix
        double *ar {nullptr}; //!< Non-zero values of the LP matrix
        double *Value {nullptr}; //!< Solution values of the LP variables
        unsigned int sizeCol {0}; //!< Number of columns (variables)
        unsigned int sizeRow {0}; //!< Number of rows (constraints)
        unsigned int sizeVet {0}; //!< Total number of non-zero elements in the matrix
        bool solved {false}; //!< Flag indicating if the LP has been solved
        string filename {""}; //!< Filename of the LP problem
        unordered_map<string, unsigned int> ReactionsNamesId; //!< Mapping of reaction names to numeric IDs
        vector<string> ReactionsNamesOrd; //!< Ordered list of reaction names
        ifstream in_var; //!< Input stream for variability data
        ofstream out_var; //!< Output stream for variability results
        unsigned int flux_var; //!< Index of the flux variable considered for variability
        double gamma {1.0}; //!< Gamma value for variability

        //! Helper function to set the type of bounds based on string input
        int setTypeBound(string typeString);

        //! Helper function to set the optimization direction (min/max) based on string input
        int setTypeObj(string typeString){
            int type;
            if (typeString.compare("GLP_MAX")==0)
                type=GLP_MAX;
            else
                type=GLP_MIN;
            return type;
        }

        //! Parses flux names from input file
        void parseFluxNames(ifstream& in, general::Parser& parser, const char* delimC);

        //! Parses model dimensions and type from input file
        void parseModelDimensionsAndType(ifstream& in, general::Parser& parser, const char* delimC, unsigned int& sizeRow, unsigned int& sizeCol, int& typeOBJ, int variability);

        //! Parses matrix size from input file
        void parseSizeVet(ifstream& in, general::Parser& parser, const char* delimC, unsigned int& sizeVet);

        //! Parses objective coefficients from input file
        void parseObjectiveCoefficients(ifstream& in, general::Parser& parser, const char* delimC, int setDefaultCoefficients, int variability, int flux_var, string* var_obj_eq);

        //! Sets row bounds from input data
        void setRowBounds(ifstream& in, general::Parser& parser, const char* delimC);

        //! Sets column bounds from input data
        void setColumnBounds(ifstream& in, general::Parser& parser, const char* delimC);

        //! Sets the sparse matrix elements from input data
        void setSparseMatrix(ifstream& in, general::Parser& parser, const char* delimC);

        //! Initializes the LP problem from file
        void initializeLP(const char* fileProb);

        //! Allocates memory for matrix indices and values
        void allocateMemory(unsigned int sizeVet);

        //! Handles variability management based on input data
        void manageVariability(ifstream& in, general::Parser& parser, const char* delimC, int variability, const char* FluxName);

        //! Adds a row for handling variability in the objective
        void addVariabilityRow(const string& var_obj_eq, general::Parser& parser, const char* delimC, unsigned int& sizeVet);

    public:
        LPprob(); //!< Default constructor
         //!< Copy constructor
        LPprob(const LPprob& t){
            if (t.filename!=""){
                this->updateLP(t.filename.c_str());
            }
        };
        LPprob(const char * FileProb); //!< Constructor taking problem file
        LPprob(const char * FileProb, const char* FileInVar, const char* FileOutVar, int typeOBJ, const char* FluxName, int gamma); //!< Full constructor for variability handling

        //! Updates the LP problem from a file
        void updateLP(const char * FileProb, int variability = 0, int typeOBJ = -1, const char* FluxName = "");

        //! Solves the LP problem using simplex method
        void solve(){
            cout<<"\n\n-------------------------------------------------------"<<endl;
            glp_simplex(lp, NULL);
            solved=true;
            cout<<"-------------------------------------------------------\n"<<endl;
        };
        //! Handles solution variability
        void solveVariability();

        //! Returns the objective function value
        inline double getOBJ(){
            if (!solved) solve();
            return glp_get_obj_val(lp);
        };

        //! Returns the solution values for variables
        inline double* getVariables(){
            if (!solved) solve();
            for (unsigned int i=1;i<=sizeCol;++i){
            Value[i]= glp_get_col_prim(lp, i);
            }
            return Value;
        };

        //! Returns the lower bound for a specified variable
        inline double getLwBounds(int indexR){
						double LB = glp_get_col_lb(lp, indexR);
            return LB;
        };

        //! Returns the upper bound for a specified variable
        inline double getUpBounds(int indexR){
	    			double UB = glp_get_col_ub(lp, indexR);
            return UB;
        };

        //! Prints the last GLPK solution to standard output
        void print(){
            if (!solved) solve();
            cout<<"Obj value:"<< getOBJ()<<endl<<endl;
            getVariables();
           auto it=ReactionsNamesOrd.begin();
            for (unsigned int i=1;i<=sizeCol;++i,++it){
                cout<<*it<<":"<<Value[i]<<endl;

            }
        };

        //! Writes variable values to an output file stream
        inline void printValue(ofstream& out){
            getVariables();
            for (unsigned int i=1;i<=sizeCol;++i){
                out<<" "<<Value[i];
            }
        };

        //! Writes the objective value to an output file stream
        inline void printObject(ofstream& out){
            out<<" "<<getOBJ();
        }

        //! Writes the upper bounds of variables to an output file stream
        inline void printUpper(ofstream& out){
            for (unsigned int i=1;i<=sizeCol;++i){
                out<<" "<<glp_get_col_ub(lp, i);
            }

        }

        //! Writes the lower bounds of variables to an output file stream
        inline void printLower(ofstream& out){
            for (unsigned int i=1;i<=sizeCol;++i){
                out<<" "<<glp_get_col_lb(lp, i);
            }

        }

        //! Writes the bounds of variables to an output file stream
        inline void printLowerMax(ofstream& out){
            for (unsigned int i=1;i<=sizeCol;++i){
                out<<" "<<glp_get_col_lb(lp, i)<<" "<<glp_get_col_ub(lp, i);
            }

        }

        //! Writes the flux names to an output file stream
        inline void printFluxName(ofstream& out){
            for (auto it=ReactionsNamesOrd.begin();it!=ReactionsNamesOrd.end();++it){
                out<<" "<<*it;
            }
        };

        //! Writes flux names with "_Lb" and "_Ub" suffixes to an output file stream
        inline void printFluxNameMinMax(ofstream& out){
            for (auto it=ReactionsNamesOrd.begin();it!=ReactionsNamesOrd.end();++it){
                out<<" "<<*it<<"_Lb"<<" "<<*it<<"_Ub";
            }
        };


        //! Updates the bounds for a specified variable
        inline void update_bound(int indexR, string TypeBound, double Lb, double Ub){
            glp_set_col_bnds(lp, indexR, setTypeBound(TypeBound) , Lb, Ub);
            cout<<"Bounds of "<< indexR <<" is updated as: ["<<Lb<<";"<<Ub<<"]"<<endl;
        };

        inline void update_bound(int indexR, int TypeBound, double Lb, double Ub){
            glp_set_col_bnds(lp, indexR, TypeBound , Lb, Ub);          
            //cout<<"Bounds of "<< indexR <<" is updated as: ["<<Lb<<";"<<Ub<<"]"<<endl;         
        };

        //!Returns the type of j-th column, i.e. the type of corresponding structural variable, as follows: GLP_FR — free (unbounded) variable; GLP_LO — variable with lower bound; GLP_UP — variable with upper bound; GLP_DB — double-bounded variable; GLP_FX — fixed variable.
        inline int  get_bound_type(int indexR){
            return glp_get_col_type(lp, indexR);
        }

        //! Converts a flux name to its numeric ID
        inline int fromNametoid(const string& name){
            auto it=ReactionsNamesId.find(name);
            if (it!=ReactionsNamesId.end())
                return it->second;
            else
                return -1;

        }

        //! Writes the solution to a specified file
				inline void writeSolutionToFile(const std::string& filename) {
						if (!solved) {
								throw Exception("LP problem not solved yet.");
						}

						std::ofstream outFile(filename);
						if (!outFile) {
								throw Exception("Could not open file for writing: " + filename);
						}

						// Write the objective value
						outFile << "Objective Value: " << getOBJ() << std::endl;

						// Fetch and write each variable's value
						getVariables();  // Ensure variables are updated
						outFile << "Variable Values:" << std::endl;
						for (unsigned int i = 1; i <= sizeCol; ++i) {
								outFile << ReactionsNamesOrd[i - 1] << ": " << Value[i] << std::endl;
						}

						outFile.close();
				}

        ~LPprob() {
            //--count;
            if (sizeVet){
            free(ia);
            free(ja);
            free(ar);
            glp_delete_prob(lp);
            //if (count==0)
            //    glp_free_env();
            };
        };
    };
}


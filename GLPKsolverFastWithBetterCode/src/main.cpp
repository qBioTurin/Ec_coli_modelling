#include "GLPKsolve.hpp"   // Include the header for the GLPK solver class definitions.
#include <iostream>        // Include for input and output stream objects, e.g., std::cout.
#include <sstream>         // Include for string stream objects (std::stringstream).
#include <string>          // Include for std::string class.
#include <chrono>          // Include for high-resolution clock and duration classes.
#include <sys/resource.h>  // Include for retrieving resource usage (memory, CPU time).
#include <functional>      // Include for std::function, used to pass functions as arguments.

using namespace FBGLPK;   // Use the FBGLPK namespace to avoid prefixing its identifiers.
using namespace std::chrono; // Use the std::chrono namespace for easy access to time-related classes.

// Function to measure the execution time of a given function.
milliseconds measureExecutionTime(std::function<void()> func) {
    auto start = high_resolution_clock::now();  // Record start time
    func();  // Execute the function
    auto stop = high_resolution_clock::now();   // Record stop time
    return duration_cast<milliseconds>(stop - start);  // Return the duration in milliseconds
}

// Function to measure the maximum memory usage of the current process.
long measureMemoryUsage() {
    struct rusage usage;  // Declare a structure to hold resource usage statistics
    getrusage(RUSAGE_SELF, &usage);  // Get resource usage for the calling process
    return usage.ru_maxrss;  // Return max resident set size used (in kilobytes)
}

// Function to run the solver and capture its console output, time, and memory usage.
std::string runSolver(const std::string& inputFilePath) {
    std::stringstream capture;  // Create a stringstream to capture output
    std::streambuf* old = std::cout.rdbuf(capture.rdbuf());  // Redirect std::cout to our stringstream

    LPprob lpProblem(inputFilePath.c_str());  // Create an LP problem instance from the file path
    auto duration = measureExecutionTime([&]() { lpProblem.solve(); });  // Measure time to solve LP
    long memoryUsed = measureMemoryUsage();  // Measure memory usage after solving

    // Reset std::cout to its old self
    std::cout.rdbuf(old);

    // Prepare the final output string with results and statistics
    std::string output = capture.str();  // Get all captured output
    output += "Time taken to solve LP: " + std::to_string(duration.count()) + " milliseconds\n";
    output += "Memory used: " + std::to_string(memoryUsed) + " kilobytes\n";
    return output;
}

// Main function: Entry point of the program
int main(int argc, char* argv[]) {
    if (argc != 2) {  // Ensure there is exactly one command line argument (excluding program name)
        std::cerr << "Usage: " << argv[0] << " <path_to_input_file>" << std::endl;
        return 1;  // Exit with an error code
    }

    std::string results = runSolver(argv[1]);  // Run the solver using the provided file path
    std::cout << results;  // Print the results to the terminal
    return 0;  // Exit successfully
}


// STL includes
#include <iostream>
#include <fstream>
#include <set>
#include <cstdlib> // error codes
#include <memory>
#include <cmath>
#include <chrono>
#include <random>
// Local includes
#include "models.h"
#include "types.h"
#include "config.h"

#if HAVE_LIBBOOST_PROGRAM_OPTIONS == 1
  // Boost
  #include <boost/program_options.hpp>
  namespace po = boost::program_options;
#endif


std::string model_desc_help_string(std::set<std::string> models,
                                   std::map<std::string, std::string> models_desc,
                                   std::string base_string)
{
  std::string help_string = base_string;
  for (auto m: models)
    help_string += " \"" + m + "\"" + " : " + models_desc.at(m) + "\n";
  return help_string.substr(0, help_string.size() - 1);
}


int main(int argc, char const *argv[])
{
  // models
  std::map<std::string, std::string> models_description;
  std::map<std::string, std::string> models_parameters;
  std::map<std::string, unsigned int> models_parameters_count;
  std::set<std::string> implemented_models = {"delayed", "gn", "bianconi_barabasi", "generalized_gn"};
  // unstructured models
  models_description["delayed"] = "L. H.-Dufresne. et al., PRE 92, (2015).";
  models_parameters["delayed"] = "a, alpha, b, mu, tau.";
  models_parameters_count["delayed"] = 5;
  // structured models
  models_description["gn"] = "Krapivsky-Redner, PRE 63, (2001).";
  models_parameters["gn"] = "gamma, m.";
  models_parameters_count["gn"] = 2;
  models_description["bianconi_barabasi"] = "Bianconie-Barabasi, PRL 86, (2001).";
  models_parameters["bianconi_barabasi"] = "gamma, nu, avg, std, m.";
  models_parameters_count["bianconi_barabasi"] = 5;
  models_description["generalized_gn"] = "General PA with densification.";
  models_parameters["generalized_gn"] = "a, alpha, b, m1, m2, gamma, tau, negative_kernel.";
  models_parameters_count["generalized_gn"] = 8;

  /* ~~~~~ Program options ~~~~~~~*/
  unsigned int T;
  double_vec_t parameters;
  std::string model_name;
  unsigned int seed = 0;
  std::string state_file_path;



  // ================
  // BOOST INTERFACE
  // ================
  #if HAVE_LIBBOOST_PROGRAM_OPTIONS == 1
    // boost::po call
    po::options_description description("Options");
    description.add_options()
    ("model,m", po::value<std::string>(&model_name)->default_value("gn"),
      // automatically generated help message
      model_desc_help_string(implemented_models,
                             models_description,
                             "Name of the growth model. The implemented models are:\n").c_str())
    ("parameters,p", po::value<double_vec_t>(&parameters)->multitoken(),
      // automatically generated help message
      model_desc_help_string(implemented_models,
                             models_parameters,
                             "Parameters of the model, provided as a list of doubles.\n").c_str())
    ("T,t", po::value<unsigned int>(&T)->default_value(1000), "Number of growth events.")
    ("seed,d", po::value<unsigned int>(&seed),
        "Seed of the pseudo random number generator (Mersenne-twister 19937)."\
        "Seeded with current time if seed is not specified or equal to 0.")
    ("state_file,f", po::value<std::string>(&state_file_path),
      "Save state variables in a file upon completion (degree, fitness, etc.).")
    ("verbose,v", "Output parameters to stdlog.")
    ("help,h", "Produce this help message.")
    ;
    po::variables_map var_map;
    try
    {
      po::store(po::parse_command_line(argc,argv,description), var_map);
      po::notify(var_map);
    }
    catch (po::validation_error& e)
    {
      std::clog << "Boost program option error:\n";
      std::clog << e.what();
      std::clog << "\n";
      return EXIT_FAILURE;
    }

    // Input validation and actions
    if (var_map.count("help") > 0 || argc == 1)
    {
      std::clog << "Usage:\n"
                << "  "+std::string(argv[0])+" [--option_1=value] [--option_s2=value] ...\n";
      std::clog << description;
      return EXIT_SUCCESS;
    }
    if (var_map.count("seed") == 0)
    {
      seed = (unsigned int) std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    if (models_parameters_count[model_name] != parameters.size())
    {
      std::clog << "Incorrect number of parameters for the \"" + model_name + "\" growth model.\n";
      return EXIT_FAILURE;
    }
    if (implemented_models.find(model_name) == implemented_models.end())
    {
      std::clog << "Model \"" + model_name + "\" not implemented.\n";
      return EXIT_FAILURE;
    }
    if (T <= 1)
    {
      std::clog << "Number of event T is too small (T=" << T << ")\n";
      return EXIT_FAILURE;
    }

    // Logger
    if (var_map.count("verbose") > 0)
    {
      std::clog << "Model: " << model_name << "\n";
      std::clog << "Parameters (" + models_parameters[model_name] + "): ";
      for (auto p: parameters)
        std::clog << p << " ";
      std::clog << "\n";
      std::clog << "T: " << T << "\n";
      std::clog << "Seed: " << seed << "\n";
    }
  // ================
  // end of boost interface
  // ================
  #else
    model_name = std::string(argv[1]);
    if (model_name=="-h")
    {
      std::clog << "This is the limited interface of this program. boost::program_options could not be found and linked.\n";
      std::clog << "Usage: " << argv[0] << "  model_name T seed param1 param2 ...\n\n"; 
      std::clog << model_desc_help_string(implemented_models,
                                          models_description,
                                          "Name of the growth model. The implemented models are:\n");
      std::clog << "\n\n";
      std::clog << model_desc_help_string(implemented_models,
                                          models_parameters,
                                          "Parameters of the model, provided as a list of doubles.\n");
      std::clog << "\n";
      return EXIT_SUCCESS;
    }
    T = std::atoi(argv[2]);
    seed = std::atoi(argv[3]);
    parameters.resize(models_parameters_count[model_name], 0);
    for (unsigned int i = 0; i < models_parameters_count[model_name]; ++i)
    {
      parameters[i] = std::atof(argv[4 + i]);
    }
  #endif

  /*~~~~~~~Initialize model~~~~~~~~~~~~*/
  std::shared_ptr<growth_model> model;

  // unstructured models
  if (model_name == "delayed")
  {
    double a = parameters[0];
    double alpha = parameters[1];
    double b = parameters[2];
    double gamma = parameters[3];
    double tau = parameters[4];
    model = std::make_shared<delayed_model>(a, alpha, b, gamma, tau);
  }
  // structured models
  if (model_name == "gn")
  {
    double gamma = parameters[0];
    unsigned int m = (unsigned int) parameters[1];
    model = std::make_shared<gn_model>(gamma, m);
  }
  if (model_name == "generalized_gn")
  {
    double a = parameters[0];
    double alpha = parameters[1];
    double b = parameters[2];
    unsigned int m1 = (unsigned int) parameters[3];
    unsigned int m2 = (unsigned int) parameters[4];
    double gamma = parameters[5];
    double tau = parameters[6];
    if (parameters[7]) gamma = -gamma;
    model = std::make_shared<generalized_gn_model>(a, alpha, b, m1, m2, gamma, tau);
    --T; // we start with an edge
  }
  if (model_name == "bianconi_barabasi")
  {
    double gamma = parameters[0];
    double nu = parameters[1];
    double avg = parameters[2];
    double std = parameters[3];
    unsigned int m = (unsigned int) parameters[4];
    model = std::make_shared<bianconi_barabasi_model>(gamma, nu, avg, std, m);
  }


  /*~~~~~~~~~~~Run~~~~~~~~~~~~~~~~*/
  std::mt19937 engine(seed);
  model->run(T, engine);
  // std::clog << "Ran\n";

  /*~~~Dump history to stream~~~~*/
  model->print_history(std::cout);
  #if HAVE_LIBBOOST_PROGRAM_OPTIONS == 1
      if (var_map.count("state_file") > 0)
      {
        std::ofstream sf(state_file_path.c_str(), std::ios::out);
        model->print_states(sf);
        sf.close();
      }
  #endif


  return 0;
}

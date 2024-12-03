#' Generate Julia Code for ODE Models
#'
#' This function generates Julia code to solve ordinary differential equation (ODE) models
#' \deqn{\displaystyle \frac{dx}{dt} = f(x, p)\,, \quad x(t = t_0) = x_0(p)\,,}
#' along with supporting functions for integration and Jacobian computation. The function
#' automatically translates an equation list into Julia syntax and writes the resulting code
#' to a file. It also creates an interface for solving ODEs and computing the Jacobian
#' \deqn{\displaystyle \frac{\partial x}{\partial p}}
#' from R using the JuliaConnectoR package.
#'
#' @param odefunction A named list where names are the dynamic variables and values are the corresponding ODEs written as strings.
#' @param modelname A character string specifying the name of the generated Julia model function. Default is `"odemodel"`.
#' @param file A character string specifying the filename where the generated Julia code will be saved. Default is constructed as `paste0(modelname, ".jl")`.
#' @param events An optional data frame specifying events to occur at certain times. It must contain the columns: `"var"`, `"time"`, `"value"`, and `"method"`. The `"value"` column can contain expressions involving dynamic variables and parameters (e.g., `"0.5 * A"`, `"k2"`).
#'
#' @return An object with attributes:
#' \itemize{
#'   \item `"equations"`: The ODE equations passed to the function.
#'   \item `"variables"`: The dynamic variables in the system.
#'   \item `"parameters"`: The parameters in the system.
#'   \item `"events"`: The events data frame, if provided.
#'   \item `"modelname"`: The name of the generated Julia model.
#'   \item `"juliacode"`: The generated Julia code.
#' }
#'
#' The returned object contains two methods:
#' \itemize{
#'   \item `$solve(x0, dynpars, times, solver = "AutoTsit5(Rosenbrock32())", atol = 1e-8, rtol = 1e-6, maxsteps = 1e5)`: Solves the ODE system.
#'   \item `$senssolve(x0, dynpars, times, solver = "AutoTsit5(Rosenbrock32())", atol = 1e-8, rtol = 1e-6, maxsteps = 1e5))`: Solves the ODE system and computes sensitivities (Jacobian).
#' }
#'
#' @examples
#' \dontrun{
#' # Define ODE system
#' odefunction <- list(
#'   Prey = "alpha * Prey - beta * Prey * Predator",
#'   Predator = "delta * Prey * Predator - gamma * Predator"
#' )
#'
#' # Define events
#' events <- data.frame(
#'   var = c("Prey", "Predator"),
#'   time = c(50, 100),
#'   value = c("2 * Prey", "Predator / 2"),
#'   method = c("add", "rep")
#' )
#'
#' # Generate the Julia model
#' odemodel <- juliaODEmodel(odefunction, modelname = "LotkaVolterra", events = events)
#'
#' # Initial conditions and parameters
#' inits <- c(Prey = 40, Predator = 9)
#' params <- c(alpha = 0.1, beta = 0.02, delta = 0.01, gamma = 0.1)
#' times <- seq(0, 200, length.out = 500)
#'
#' # Solve the ODE system without sensitivities
#' solution <- odemodel$solve(inits, params, times)
#' print(head(solution))
#'
#' # Solve the ODE system with sensitivities
#' solution_sens <- odemodel$senssolve(inits, params, times)
#' print(head(solution_sens))
#' }
#' @export
juliaODEmodel <- function(odefunction, modelname = "odemodel", file = paste0(modelname, ".jl"), events = NULL) {

  # ================================
  # 1. Define Reserved Symbols
  # ================================
  # 't' and 'time' are treated as time variables and not as parameters
  reserved_symbols <- c("t", "time")

  # ================================
  # 2. Extract Dynamic Variables and Symbols from ODEs
  # ================================
  dynvars <- names(odefunction)
  symbols_ode <- unique(unlist(lapply(odefunction, function(expr) {
    parsed_expr <- parse(text = expr)
    all.vars(parsed_expr)
  })))

  # ================================
  # 3. Extract Symbols from Events' Value Expressions
  # ================================
  if (!is.null(events)) {
    symbols_events <- unique(unlist(lapply(events$value, function(expr) {
      parsed_expr <- parse(text = expr)
      all.vars(parsed_expr)
    })))
  } else {
    symbols_events <- character(0)
  }

  # ================================
  # 4. Determine All Parameters
  # ================================
  # Parameters are symbols that are neither dynamic variables nor reserved symbols
  all_symbols <- union(symbols_ode, symbols_events)
  parameters <- setdiff(all_symbols, c(dynvars, reserved_symbols))

  # Indices for dynamic variables and parameters
  dynvars_indices <- stats::setNames(seq_along(dynvars), dynvars)
  parameter_indices <- stats::setNames(seq_along(parameters), parameters)

  # ================================
  # 5. Helper Function to Replace Symbols in Expressions
  # ================================
  replace_symbols <- function(expr) {
    # Replace dynamic variables with u[index]
    for (x in dynvars) {
      expr <- gsub(paste0("\\b", x, "\\b"), paste0("u[", dynvars_indices[x], "]"), expr)
    }
    # Replace parameters with p[index]
    for (pname in parameters) {
      expr <- gsub(paste0("\\b", pname, "\\b"), paste0("p[", parameter_indices[pname], "]"), expr)
    }
    # Replace 'time' and 't' with 't'
    expr <- gsub("\\b(time|t)\\b", "t", expr)
    return(expr)
  }

  # ================================
  # 6. Generate Julia Code for ODE Equations
  # ================================
  eqncode <- lapply(odefunction, replace_symbols)

  # Combine the equations into the Julia ODE function
  odefuncode <- paste0(
    "# Define the evolution function \n",
    "function ", modelname, "!(du, u, p, t)\n",
    paste0("    du[", seq_along(dynvars), "] = ", eqncode, collapse = "\n"),
    "\nend \n\n"
  )

  # ================================
  # 7. Generate Julia Code for Events (if provided)
  # ================================
  if (!is.null(events)){
    # Check for required columns in events
    if (!all(c("var", "time", "value", "method") %in% colnames(events))) {
      stop("The events dataframe must contain columns: 'var', 'time', 'value', 'method'.")
    }

    # Initialize code blocks for conditions and affects
    affect_code <- NULL
    condition_code <- NULL
    cb_code <- NULL

    # Iterate over each event to generate corresponding Julia code
    for (row in seq_len(nrow(events))) {
      var <- events$var[row]
      time <- events$time[row]
      value_expr <- events$value[row]
      method <- events$method[row]

      # Define the condition for the event
      condition_code <- paste0(condition_code,
                               "condition_", as.character(row), "(u, t, integrator) = integrator.t == ", as.character(time), "\n")

      # ================================
      # 5.a. Replace Symbols in Event Value Expressions
      # ================================
      # Start with the original value expression
      value_replaced <- value_expr

      # Replace dynamic variables with integrator.u[index]
      for (x in dynvars) {
        value_replaced <- gsub(paste0("\\b", x, "\\b"), paste0("integrator.u[", dynvars_indices[x], "]"), value_replaced)
      }

      # Replace parameters with integrator.p[index]
      for (pname in parameters) {
        value_replaced <- gsub(paste0("\\b", pname, "\\b"), paste0("integrator.p[", parameter_indices[pname], "]"), value_replaced)
      }

      # Replace 'time' and 't' with integrator.t
      value_replaced <- gsub("\\b(time|t)\\b", "integrator.t", value_replaced)

      # ================================
      # 5.b. Map the Method to the Corresponding Julia Operation
      # ================================
      operation_code <- switch(
        method,
        add = paste0("integrator.u[", dynvars_indices[var], "] += ", value_replaced),
        mult = paste0("integrator.u[", dynvars_indices[var], "] *= ", value_replaced),
        rep = paste0("integrator.u[", dynvars_indices[var], "] = ", value_replaced),
        stop(paste("Unsupported method:", method))
      )

      # Define the function that affects the system state
      affect_code <- paste0(affect_code,
                            "affect_", as.character(row), "!(integrator) = ", operation_code, "\n")

      # Add the callback for this event
      cb_code <- paste0(cb_code,
                        "cb_", as.character(row), " = DiscreteCallback(condition_", as.character(row),
                        ", affect_", as.character(row), "!, save_positions = (true, true)) \n")
    }

    # Combine all callbacks into a CallbackSet
    cb_code <- paste0(cb_code,
                      "cbs = CallbackSet(", paste0("cb_", seq_along(events$time), collapse = ", "), ")")

    # Combine all event-related Julia code
    event_code <- paste0(
      "# Julia functions for event handling \n",
      condition_code, "\n",
      affect_code, "\n",
      cb_code, "\n",
      "tstop_events = ", paste0("[", paste(unique(events$time), collapse = ", "), "]"), "\n"
    )

  } else{
    # If no events are defined
    event_code <- paste0("# Julia functions for event handling \n",
                         "cbs = nothing \n",
                         "tstop_events = nothing \n\n")
  }


  # ================================
  # 8. Combine All Julia Code
  # ================================
  julia_code <- paste0(
    "## Autogenerated Julia Code by JuliaODEmodel --- \n\n",
    "# Loading necessary Julia packages \n",
    "using OrdinaryDiffEq, ForwardDiff\n\n",
    odefuncode,
    event_code, "\n\n",

    "# Function for solving the ODEs without sensitivities\n",
    "function solve_", modelname, "(p, times, solver::SciMLBase.AbstractODEAlgorithm, args...; kwargs...)\n",
    "    sol = solve(ODEProblem(", modelname, "!, (@view p[1:", as.character(length(dynvars)), "]), ",
    "(times[1], times[end]), (@view p[", as.character(length(dynvars)+1), ":end])), ",
    "solver, args...; saveat=times, kwargs...)\n",
    "    return hcat(sol.t, Array(sol)')\n",
    "end \n\n",

    "# Function for solving the ODEs with sensitivities\n",
    "function solvesens_", modelname, "(p, times, solver::SciMLBase.AbstractODEAlgorithm, args...; kwargs...)\n",
    "    sol = solve(ODEProblem(", modelname, "!, (@view p[1:", as.character(length(dynvars)), "]), ",
    "(times[1], times[end]), (@view p[", as.character(length(dynvars)+1), ":end])), ",
    "solver, args...; saveat=times, kwargs...)\n",
    "    jac = reshape(ForwardDiff.jacobian(pars -> hcat(solve(ODEProblem(", modelname, "!, (@view pars[1:",
    as.character(length(dynvars)), "]), (times[1], times[end]), (@view pars[",
    as.character(length(dynvars)+1), ":end])), solver, args...; saveat=times, kwargs...).u...)', p), ",
    "length(sol.t), :)\n",
    "    return hcat(sol.t, Array(sol)', jac)\n",
    "end \n\n"
  )

  # ================================
  # 9. Save and Execute Julia Code
  # ================================
  # Save the generated Julia code to the specified file
  writeLines(julia_code, file)

  # Execute the Julia code in the Julia session
  JuliaConnectoR::juliaEval(julia_code)

  # ================================
  # 10. Define the ODEmodel Object with Solution Methods
  # ================================
  ODEmodel <- list()

  # Function to solve ODEs without sensitivities
  ODEmodel$solve <- function(x0, dynpars, times, solver = "AutoTsit5(Rosenbrock32())", atol = 1e-8, rtol = 1e-6, maxsteps = 1e5) {
    # Integrate in Julia
    out <- JuliaConnectoR::juliaLet(
      paste0("solve_", modelname,"([x0;dynpars], times, ", solver, "; callback = cbs, tstops = tstop_events, abstol = atol, reltol = rtol, maxiters = maxsteps)"),
      x0 = x0, dynpars = dynpars, times = times, atol = atol, rtol = rtol, maxsteps = maxsteps
    )
    colnames(out) <- c("time", dynvars)
    return(out)
  }

  # Function to solve ODEs with sensitivities
  ODEmodel$senssolve <- function(x0, dynpars, times, solver = "AutoTsit5(Rosenbrock32())", atol = 1e-8, rtol = 1e-6, maxsteps = 1e5) {
    # Calculate Jacobian in Julia
    out <- JuliaConnectoR::juliaLet(
      paste0("solvesens_", modelname,"([x0;dynpars], times, ", solver, "; callback = cbs, tstops = tstop_events, abstol = atol, reltol = rtol, maxiters = maxsteps)"),
      x0 = x0, dynpars = dynpars, times = times, atol = atol, rtol = rtol, maxsteps = maxsteps
    )

    # Generate column names for the Jacobian
    param_names <- c(dynvars, parameters)
    jac_colnames <- c("time",dynvars, as.vector(outer(dynvars, param_names, paste, sep = ".")))
    colnames(out) <- jac_colnames
    return(out)
  }

  # ================================
  # 11. Add Attributes to the ODEmodel Object
  # ================================
  attr(ODEmodel, "equations") <- odefunction
  attr(ODEmodel, "variables") <- dynvars
  attr(ODEmodel, "parameters") <- parameters
  attr(ODEmodel, "events") <- events
  attr(ODEmodel, "modelname") <- modelname
  attr(ODEmodel, "juliacode") <- julia_code

  # ================================
  # 12. Return the ODEmodel Object
  # ================================
  return(ODEmodel)
}




#' .onAttach: Ensure Julia Dependencies on Package Load
#'
#' This function is called automatically when the package is attached to the R session
#' (e.g., using `library(yourpackage)`).
#'
#' @param libname The library path to the package.
#' @param pkgname The name of the package being attached.
#'
#' @details
#' The `.onAttach` function ensures that all required Julia dependencies are installed
#' by calling `ensureJuliaPackages`. It provides startup messages to inform the user
#' about the status of the installation.
#'
#' @note
#' Users do not need to call this function directly; it is executed automatically.
#'
#' @seealso
#' \code{\link{ensureJuliaPackages}}
#'
#' @export
ensureJuliaPackages <- function() {
  # List of required Julia packages
  required_packages <- c("OrdinaryDiffEq", "ForwardDiff")

  # Fetch the list of installed packages and extract names
  installed_packages <- JuliaConnectoR::juliaEval("
    [info.name for (uuid, info) in Pkg.dependencies()]
  ")

  # Check and install missing packages
  for (pkg in required_packages) {
    if (!(pkg %in% installed_packages)) {
      packageStartupMessage(sprintf("Installing missing Julia package: %s", pkg))
      JuliaConnectoR::juliaCall("Pkg.add", pkg)
    } else {
      packageStartupMessage(sprintf("Julia package %s is already installed.", pkg))
    }
  }

  packageStartupMessage("All required Julia packages are now installed.")
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Checking Julia dependencies...")
  ensureJuliaPackages() # Call the function to ensure required Julia packages are installed.
  packageStartupMessage("Julia dependencies are ready.")
}





#' Generate ODE model object that delivers Julia integration methods
#'
#' This function generates Julia code to solve ordinary differential equation (ODE) models
#' \deqn{\displaystyle \dot{x} = f(x, p)\,, \quad x(t = t_0) = x_0(p)\,,}
#' along with supporting functions for integration and Jacobian computation. The function
#' automatically translates an equation list into Julia syntax and writes the resulting code
#' to a file. It also creates an interface for solving ODEs and computing the Jacobian
#' \eqn{\partial x(t,p) /\partial p}
#' from R using the JuliaCall package.
#'
#' @param odefunction A named list where names are the dynamic variables and values are the corresponding ODEs written as strings.
#' @param modelname A character string specifying the name of the generated Julia model function. Default is "odemodel".
#' @param file A character string specifying the filename where the generated Julia code will be saved. Default is constructed as paste0(modelname, ".jl").
#' @param events An optional data frame specifying events to occur at certain times. It must contain the columns: "var", "time", "value", and "method". The "value" column can contain expressions involving dynamic variables and parameters (e.g., "0.5 * A", "k2").
#'
#' @return An object with attributes:
#' \itemize{
#'   \item "equations": The ODE equations passed to the function.
#'   \item "variables": The dynamic variables in the system.
#'   \item "sensvariables": The labels of the derivatives. E.g. Prey.alpha for \eqn{\partial \text{Prey} /\partial alpha}
#'   \item "parameters": The parameters in the system.
#'   \item "events": The events data frame, if provided.
#'   \item "modelname": The name of the generated Julia model.
#'   \item "juliacode": The generated Julia code.
#' }
#'
#' The returned object contains two methods:
#' \itemize{
#'   \item $solve(inits, dynpars, times, optionsOde = NULL): Solves the ODE system.
#'   \item $senssolve(inits, dynpars, times, optionsOde = NULL): Solves the ODE system and computes sensitivities \eqn{\partial x(t,p) /\partial p}.
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
#'   method = c("add", "replace")
#' )
#'
#' # Generate the Julia model
#' odemodel <- odemodelJL(odefunction, modelname = "LotkaVolterra", events = events)
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
odemodelJL <- function(odefunction, modelname = "odemodel", file = paste0(modelname, ".jl"), events = NULL) {

  # ================================================================
  # 1. Define Reserved Symbols
  # ================================================================
  # 't' and 'time' are treated as time variables and not as parameters
  reserved_symbols <- c("t", "time")

  # ================================================================
  # 2. Extract Dynamic Variables and Symbols from ODE
  # ================================================================
  dynvars <- names(odefunction)
  symbols_ode <- unique(unlist(lapply(odefunction, function(expr) {
    parsed_expr <- parse(text = expr)
    all.vars(parsed_expr)
  })))

  # ================================================================
  # 3. Extract Symbols from Events' Value Expressions
  # ================================================================
  if (!is.null(events)) {
    symbols_events <- unique(unlist(lapply(events$value, function(expr) {
      parsed_expr <- parse(text = expr)
      all.vars(parsed_expr)
    })))
  } else {
    symbols_events <- character(0)
  }

  # ================================================================
  # 4. Determine All Parameters
  # ================================================================
  # Parameters are symbols that are neither dynamic variables nor reserved symbols
  all_symbols <- union(symbols_ode, symbols_events)
  parameters <- setdiff(all_symbols, c(dynvars, reserved_symbols))

  # Indices for dynamic variables and parameters
  dynvars_indices <- stats::setNames(seq_along(dynvars), dynvars)
  parameter_indices <- stats::setNames(seq_along(parameters), parameters)

  # ================================================================
  # 6. Generate Julia Code for ODE Function
  # ================================================================
  eqncode <- lapply(odefunction, function(expr) {
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
  })

  # Combine the equations into the Julia ODE function
  odefuncode <- paste0(
    "# Define the evolution function \n",
    "function ", modelname, "!(du, u, p, t)\n",
    paste0("    du[", seq_along(dynvars), "] = ", eqncode, collapse = "\n"),
    "\nend \n\n"
  )

  # ================================================================
  # 7. Generate Julia Code for Events (if provided)
  # ================================================================
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

      # ================================================================
      # 5.a. Replace Symbols in Event Value Expressions
      # ================================================================
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

      # ================================================================
      # 5.b. Map the Method to the Corresponding Julia Operation
      # ================================================================
      operation_code <- switch(
        method,
        add = paste0("integrator.u[", dynvars_indices[var], "] += ", value_replaced),
        multiply = paste0("integrator.u[", dynvars_indices[var], "] *= ", value_replaced),
        replace = paste0("integrator.u[", dynvars_indices[var], "] = ", value_replaced),
        stop(paste("Unsupported event method:", method))
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


  # ================================================================
  # 8. Combine All Julia Code
  # ================================================================
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

  # ================================================================
  # 9. Save and Execute Julia Code
  # ================================================================
  writeLines(julia_code, file)
  JuliaCall::julia_source(file)

  # ================================================================
  # 10. Define the ODEmodel Object with Solution Methods
  # ================================================================
  ODEmodel <- list()

  # Function to solve ODEs without sensitivities
  ODEmodel$solve <- function(inits, dynpars, times, optionsOde = NULL) {

    options <- list(method = "AutoTsit5(Rosenbrock32())", atol = 1e-6, rtol = 1e-6, maxsteps = 1e7) # Default values
    if(!is.null(options)) options[match(names(optionsOde), names(options))] <- optionsOde

    # Ensure that inits and dynpars are in the correct order
    inits <- inits[dynvars]
    dynpars <- dynpars[parameters]

    # Assign arguments in Julia
    JuliaCall::julia_assign("inits", inits)
    JuliaCall::julia_assign("dynpars", dynpars)
    JuliaCall::julia_assign("times", times)
    JuliaCall::julia_assign("atol", options$atol)
    JuliaCall::julia_assign("rtol", options$rtol)
    JuliaCall::julia_assign("maxsteps", options$maxsteps)

    # Solve
    out <- JuliaCall::julia_eval(
      paste0("solve_", modelname, "([inits; dynpars], times, ", options$method, "; callback = cbs, tstops = tstop_events, abstol = atol, reltol = rtol, maxiters = maxsteps)")
    )

    # Return
    colnames(out) <- c("time", dynvars)
    return(out)
  }

  # Function to solve ODEs with sensitivities
  ODEmodel$senssolve <- function(inits, dynpars, times, optionsOde = NULL) {

    options <- list(method = "AutoTsit5(Rosenbrock32())", atol = 1e-6, rtol = 1e-6, maxsteps = 1e7) # Default values
    if(!is.null(options)) options[match(names(optionsOde), names(options))] <- optionsOde

    # Ensure that inits and dynpars are in the correct order
    inits <- inits[dynvars]
    dynpars <- dynpars[parameters]

    # Assign arguments in Julia
    JuliaCall::julia_assign("inits", inits)
    JuliaCall::julia_assign("dynpars", dynpars)
    JuliaCall::julia_assign("times", times)
    JuliaCall::julia_assign("atol", options$atol)
    JuliaCall::julia_assign("rtol", options$rtol)
    JuliaCall::julia_assign("maxsteps", options$maxsteps)

    # Integration and Jacobian calculation in Julia
    out_total <- JuliaCall::julia_eval(
      paste0("solvesens_", modelname,"([inits;dynpars], times, ", options$method, "; callback = cbs, tstops = tstop_events, abstol = atol, reltol = rtol, maxiters = maxsteps)")
    )

    # Generate named output elements
    colnames(out_total) <- c("time",dynvars, as.vector(outer(dynvars, c(dynvars, parameters), paste, sep = ".")))
    out <- out_total[, c("time",dynvars)]
    out_sens <- out_total[, c("time", as.vector(outer(dynvars, c(dynvars, parameters), paste, sep = ".")))]

    # Return both components as a list
    return(list(out = out, out_sens = out_sens))
  }

  # ================================================================
  # 11. Return the ODEmodel Object
  # ================================================================
  attr(ODEmodel, "equations") <- odefunction
  attr(ODEmodel, "variables") <- dynvars
  attr(ODEmodel, "sensvariables") <- as.vector(outer(dynvars, c(dynvars, parameters), paste, sep = "."))
  attr(ODEmodel, "parameters") <- parameters
  attr(ODEmodel, "events") <- events
  attr(ODEmodel, "modelname") <- modelname
  attr(ODEmodel, "juliacode") <- julia_code

  return(ODEmodel)
}

#' .onAttach: Ensure Julia Dependencies on Package Load
#'
#' This function is automatically called when the package is attached to the R session.
#' It ensures that all required Julia dependencies are installed and ready for use.
#'
#' @param libname The library name (automatically provided by R when the package is loaded).
#' @param pkgname The package name (automatically provided by R when the package is loaded).
#'
#' @details
#' The `.onAttach` function ensures that all required Julia dependencies are installed
#' by calling `ensureJuliaPackages`. It provides startup messages to inform the user
#' about the status of the installation process.
#'
#' @note
#' Users do not need to call this function directly; it is executed automatically when
#' the package is loaded using `library()` or `require()`.
#'
#' @seealso
#' \code{\link{ensureJuliaPackages}}
#'
#' @export
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Checking Julia dependencies...")
  ensureJuliaPackages() # Ensure required Julia packages are installed.
  packageStartupMessage("Julia dependencies are ready.")
}

#' Ensure Required Julia Packages Are Installed
#'
#' This function checks for and installs any required Julia packages needed by the R package.
#'
#' @details
#' The function sets up the Julia environment using the `JuliaCall` package, writes a temporary
#' Julia script with functions to manage package installation, and ensures that required Julia
#' packages are available. If any package is missing, it will be installed automatically.
#'
#' @note
#' This function is called internally by `.onAttach` and does not need to be called directly
#' by the user.
#'
#' @seealso
#' \code{\link{.onAttach}}
#'
#' @import JuliaCall
#'
#' @export
ensureJuliaPackages <- function() {
  # Set up Julia environment
  JuliaCall::julia_setup()

  # Create a temporary Julia script file
  temp_julia_script <- tempfile(fileext = ".jl")

  # Ensure the temporary file is deleted when the function exits, even if an error occurs
  on.exit({
    if (file.exists(temp_julia_script)) {
      file.remove(temp_julia_script)
    }
  }, add = TRUE)

  # Write the Julia code to the temporary file
  writeLines(
    paste0(
      "# \u00A9 Rafael Arutjunjan\n",
      "using Pkg\n",
      "AddPkg(name::AbstractString) = name in keys(Pkg.project().dependencies) || Pkg.add(name)\n\n",
      "macro secureload(Ex)\n",
      "    @assert Ex.head === :using\n",
      "    for ex in Ex.args\n",
      "        AddPkg(string(ex.args[1]))\n",
      "    end\n",
      "    Ex\n",
      "end"
    ),
    temp_julia_script
  )

  # Source the Julia script
  JuliaCall::julia_source(temp_julia_script)

  # List of required Julia packages
  required_packages <- c("OrdinaryDiffEq", "ForwardDiff")

  # Check and install missing packages
  for (pkg in required_packages) {
    packageStartupMessage(sprintf("Checking for Julia package: %s", pkg))
    JuliaCall::julia_eval(sprintf("using Pkg; @secureload using %s", pkg))
  }
}







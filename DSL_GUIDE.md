# 🏴‍☠️ Coalescent.jl DSL Guide

Ahoy there, matey! Welcome to the comprehensive guide for the new Julia DSL (Domain-Specific Language) for Coalescent.jl. This DSL provides a more idiomatic and powerful way to define coalescent simulation models directly in Julia, replacing the need for YAML configuration files.

## 🌊 Why Use the DSL?

The DSL offers several advantages over YAML-based model definition:

- **Type Safety**: Julia's type system catches errors at definition time
- **Native Julia Features**: Use functions, variables, and expressions directly
- **Better IDE Support**: Syntax highlighting, autocompletion, and error checking
- **Metaprogramming Power**: Leverage Julia's macro system for advanced features
- **Performance**: No runtime YAML parsing overhead
- **Composability**: Easily combine and modify models programmatically

## 🗺️ Quick Navigation

- [Basic Usage](#basic-usage)
- [DSL Components](#dsl-components)
- [Complete Examples](#complete-examples)
- [Comparison with YAML](#comparison-with-yaml)
- [Advanced Features](#advanced-features)
- [Simulation Usage](#simulation-usage)

## 🚢 Basic Usage

The core of the DSL is the `@model` macro, which creates a `CoalescentModel` object:

```julia
using Coalescent

# Define a simple SIR model
sir_model = @model "SIR_Example" begin
    # Parameters (the treasure in our hold!)
    @parameter β 3.0    # Transmission rate
    @parameter γ 2.0    # Recovery rate
    
    # Demes (compartments we can sample from)
    @deme I 1.0         # Infected individuals
    
    # Non-demes (important variables but not samplingable)
    @nondeme S 1e5 ode=(-β*S*I/N)    # Susceptible population
    @nondeme R 0.0 ode=(γ*I)         # Recovered population
    
    # Helper variables (computed from others)
    @helper N (S + I + R)            # Total population
    
    # Birth reactions (new infections)
    @birth I => I rate=(β*S*I/N)     # Transmission
    
    # Death reactions (recovery)
    @death I rate=(γ*I)              # Recovery from infection
    
    # Time span for simulation
    @timespan 1.0 35.0
end
```

## ⚓ DSL Components

### Parameters: `@parameter`
Define model constants and rates:

```julia
@parameter β 3.0              # Simple numeric value
@parameter γ 1.0/7.0          # Expression evaluation
@parameter complex_param log(2)/10.5  # Any Julia expression
```

### Demes: `@deme`
Define population compartments that can be sampled:

```julia
@deme I 1.0                   # Simple initial value
@deme S 1000.0 ode=(-β*S*I/N) # With explicit ODE (optional)
```

### Non-demes: `@nondeme`
Define state variables that cannot be sampled but are part of the dynamics:

```julia
@nondeme R 0.0 ode=(γ*I)      # ODE is required for non-demes
```

### Birth Reactions: `@birth`
Define how new lineages are created (forward-time births become reverse-time coalescences):

```julia
@birth I => I rate=(β*S*I/N)  # Within-compartment transmission
@birth I => E rate=(β*S*I/N)  # Between-compartment transmission
```

### Death Reactions: `@death`
Define how lineages are removed (forward-time deaths become reverse-time "births"):

```julia
@death I rate=(γ*I)           # Recovery/removal from compartment
```

### Migration Reactions: `@migration`
Define movement between compartments:

```julia
@migration E => I rate=(σ*E)  # Progression from exposed to infected
@migration I1 => I2 rate=(m*I1) # Spatial migration between patches
```

### Helper Variables: `@helper`
Define computed variables to simplify expressions:

```julia
@helper N (S + I + R)         # Total population
@helper prevalence (I/N)      # Current prevalence
```

### Time Span: `@timespan`
Set the simulation time range:

```julia
@timespan 0.0 100.0          # From time 0 to time 100
```

## 🏝️ Complete Examples

### SEIR Model (with Incubation Period)

```julia
seir_model = @model "SEIR" begin
    @parameter β 3.0      # Transmission rate
    @parameter γ1 2.0     # Incubation rate (E→I)
    @parameter γ2 2.0     # Recovery rate (I→R)
    
    @deme I 1.0           # Infected (symptomatic)
    @deme E 0.0           # Exposed (incubating)
    
    @nondeme S 1e5 ode=(-β*S*I/N)    # Susceptible
    @nondeme R 0.0 ode=(γ2*I)        # Recovered
    
    @helper N (S + E + I + R)
    
    @birth I => E rate=(β*S*I/N)     # Transmission
    @migration E => I rate=(γ1*E)    # Disease progression
    @death I rate=(γ2*I)             # Recovery
    
    @timespan 1.0 35.0
end
```

### Multi-host Model with Reservoir

```julia
reservoir_model = @model "SIR_Reservoir" begin
    # Human population parameters
    @parameter β 1.5/7.0
    @parameter γ 1.0/7.0
    @parameter N_human 1e3
    
    # Reservoir population parameters
    @parameter β_res 1.5/7.0
    @parameter N_res 1e6
    
    # Cross-species transmission
    @parameter μ 1e-4*(1.0/7.0)    # Spillover rate
    
    # Demes (both species)
    @deme I_human 0.0
    @deme I_res 1.0
    
    # Non-demes
    @nondeme S_human 1e4 ode=(-β*I_human*S_human/N_human - μ*I_res*S_human/N_human)
    @nondeme S_res 1e6 ode=(-β_res*I_res*S_res/N_res)
    
    # Within-species transmission
    @birth I_human => I_human rate=(β*I_human*S_human/N_human)
    @birth I_res => I_res rate=(β_res*I_res*S_res/N_res)
    
    # Cross-species spillover
    @birth I_res => I_human rate=(μ*I_res*S_human/N_human)
    
    # Recovery
    @death I_human rate=(γ*I_human)
    @death I_res rate=(γ*I_res)
    
    @timespan 0.0 220.0
end
```

### Metapopulation Model

```julia
metapop_model = @model "Metapopulation" begin
    # Disease parameters (same for both patches)
    @parameter β 2.0
    @parameter γ 1.0
    
    # Migration parameters
    @parameter m12 0.1    # Migration rate 1→2
    @parameter m21 0.05   # Migration rate 2→1
    
    # Patch 1
    @deme I1 5.0
    @nondeme S1 1000.0 ode=(-β*S1*I1/N1)
    @nondeme R1 0.0 ode=(γ*I1)
    
    # Patch 2
    @deme I2 1.0
    @nondeme S2 2000.0 ode=(-β*S2*I2/N2)
    @nondeme R2 0.0 ode=(γ*I2)
    
    # Helpers for patch populations
    @helper N1 (S1 + I1 + R1)
    @helper N2 (S2 + I2 + R2)
    
    # Within-patch transmission
    @birth I1 => I1 rate=(β*S1*I1/N1)
    @birth I2 => I2 rate=(β*S2*I2/N2)
    
    # Between-patch migration
    @migration I1 => I2 rate=(m12*I1)
    @migration I2 => I1 rate=(m21*I2)
    
    # Recovery
    @death I1 rate=(γ*I1)
    @death I2 rate=(γ*I2)
    
    @timespan 0.0 100.0
end
```

## 📊 Comparison with YAML

### YAML Version:
```yaml
modelname: SIR
parameters:
  - name: beta
    value: 3.0
  - name: gamma
    value: 2.0

dynamic_variables:
  - name: I
    initial_value: 1.0
  - name: S
    initial_value: 1e5
    ode: -beta*S*I/N
  - name: R
    initial_value: 0.0
    ode: gamma*I

helpers:
  - name: N
    definition: S + I + R

births:
  - source: I
    recipient: I
    rate: beta*S*I/N

deaths:
  - deme: I
    rate: gamma*I

time:
  initial: 1.0
  final: 35.0
```

### DSL Version:
```julia
sir_model = @model "SIR" begin
    @parameter β 3.0
    @parameter γ 2.0
    
    @deme I 1.0
    @nondeme S 1e5 ode=(-β*S*I/N)
    @nondeme R 0.0 ode=(γ*I)
    
    @helper N (S + I + R)
    
    @birth I => I rate=(β*S*I/N)
    @death I rate=(γ*I)
    
    @timespan 1.0 35.0
end
```

**Advantages of DSL:**
- More concise and readable
- Type safety and error checking
- IDE support with syntax highlighting
- Native Julia expressions
- Easier to compose and modify programmatically

## 🔧 Advanced Features

### Programmatic Model Construction

```julia
# Create models programmatically
function create_sir_variant(β_val, γ_val, pop_size)
    return @model "SIR_Variant" begin
        @parameter β $β_val
        @parameter γ $γ_val
        
        @deme I 1.0
        @nondeme S $pop_size ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ*I)
        
        @helper N (S + I + R)
        
        @birth I => I rate=(β*S*I/N)
        @death I rate=(γ*I)
        
        @timespan 0.0 50.0
    end
end

# Create multiple variants
models = [create_sir_variant(β, 2.0, 1000) for β in [1.0, 2.0, 3.0]]
```

### Complex Expressions

```julia
complex_model = @model "Complex" begin
    @parameter R0 2.5
    @parameter infectious_period 7.0
    @parameter latent_period 3.0
    
    # Derived parameters
    @parameter β (R0 / infectious_period)
    @parameter σ (1.0 / latent_period)
    @parameter γ (1.0 / infectious_period)
    
    # Using functions in expressions
    @helper seasonal_forcing (1 + 0.1 * cos(2π * t / 365.25))
    
    @deme I 1.0
    @deme E 0.0
    @nondeme S 1e5 ode=(-β*seasonal_forcing*S*I/N)
    @nondeme R 0.0 ode=(γ*I)
    
    @helper N (S + E + I + R)
    
    @birth I => E rate=(β*seasonal_forcing*S*I/N)
    @migration E => I rate=(σ*E)
    @death I rate=(γ*I)
    
    @timespan 0.0 365.25
end
```

## 🚢 Simulation Usage

### Basic Simulation

```julia
# Create model
sir_model = @model "SIR" begin
    # ... model definition ...
end

# Method 1: Using SampleConfiguration
sample_config = SampleConfiguration(confstr = \"\"\"
sample:
  - deme: I
    time: 15.0
    size: 100
\"\"\")

tree = SimTree(sir_model, sample_config)

# Method 2: Direct arrays
sample_times = fill(15.0, 100)
sample_states = fill("I", 100)

tree = SimTree(sir_model, sample_times, sample_states)
```

### Converting to ModelFGY

If you need to work with the existing ModelFGY infrastructure:

```julia
# Convert DSL model to ModelFGY
fgy_model = to_modelfgy(sir_model)

# Use with existing functions
solution = solveodes(fgy_model)
tree = SimTree(fgy_model, sample_config)
```

## 🛠️ Implementation Details

The DSL is implemented using Julia macros and consists of:

1. **`CoalescentModel`**: Core data structure holding model components
2. **`@model`**: Main macro for model definition  
3. **Component macros**: `@parameter`, `@deme`, `@nondeme`, etc.
4. **`to_modelfgy()`**: Conversion function to existing ModelFGY format
5. **Extended `SimTree` constructors**: Direct simulation from DSL models

## 🏴‍☠️ Error Handling

The DSL includes comprehensive validation:

```julia
# This will throw an error - no demes defined
@test_throws ErrorException @model "Bad_Model" begin
    @parameter β 3.0
    @nondeme S 1000.0 ode=(-β*S)
    @timespan 0.0 10.0
end

# This will throw an error - non-deme without ODE
@test_throws ErrorException @model "Bad_Model2" begin
    @parameter β 3.0
    @deme I 1.0
    @nondeme S 1000.0  # Missing ODE!
    @timespan 0.0 10.0
end
```

## 🌟 Summary

The Coalescent.jl DSL provides a powerful, type-safe, and idiomatic way to define coalescent simulation models in pure Julia. It replaces YAML-based configurations with a more flexible and maintainable approach while maintaining full compatibility with the existing simulation infrastructure.

Key benefits:
- ✅ Type safety and compile-time error checking
- ✅ Native Julia syntax and features  
- ✅ Better IDE support and tooling
- ✅ Programmatic model construction
- ✅ Full compatibility with existing Coalescent.jl infrastructure
- ✅ More concise and readable than YAML
- ✅ Leverages Julia's metaprogramming capabilities

Set sail with the new DSL and enjoy a smoother voyage through coalescent simulation modeling!

---

*"The code be more guidelines than actual rules. Welcome aboard!"* - Captain Codehook 🏴‍☠️
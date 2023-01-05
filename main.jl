using QPSReader
using Oscar

function load_mps(path)
  loaded_mps = readqps(path)
  M = Polymake.SparseMatrix{Polymake.Rational}(loaded_mps.ncon, loaded_mps.nvar)

  # get all linear forms
  mps_getfield = x -> getfield(loaded_mps, x)
  for (row, col, val) in zip(mps_getfield.([:arows, :acols, :avals])...)
      M[row, col] = val
  end

  # get all bounds on forms
  lower_bounds = Polymake.SparseVector{Polymake.Rational}(loaded_mps.ncon)
  upper_bounds = Polymake.SparseVector{Polymake.Rational}(loaded_mps.ncon)

  for (index, (u_bound, l_bound)) in enumerate(zip(mps_getfield.([:ucon, :lcon])...))
      if l_bound != -Inf
          lower_bounds[index] = -l_bound
      end

      if u_bound != Inf
          upper_bounds[index] = u_bound
      end
  end

  M_l = hcat(lower_bounds, M)
  M_u = hcat(upper_bounds, M)
  M = vcat(M_l, M_u)


  # get all bounds on variables
  U = Polymake.to_sparsematrix_rational(Polymake.common.unit_matrix(loaded_mps.nvar))

  lower_bounds = Polymake.SparseVector{Polymake.Rational}(loaded_mps.nvar)
  upper_bounds = Polymake.SparseVector{Polymake.Rational}(loaded_mps.nvar)

  for (index, (u_bound, l_bound)) in enumerate(zip(mps_getfield.([:uvar, :lvar])...))
      if l_bound != -Inf
          lower_bounds[index] = -l_bound
      end

      if u_bound != Inf
          upper_bounds[index] = u_bound
      end
  end

  U_l = hcat(lower_bounds, U)
  U_u = hcat(upper_bounds, U)
  U = vcat(U_l, U_u)

  M = vcat(M, U)

  h_description = Polyhedron(Polymake.polytope.Polytope(INEQUALITIES=M))
  objective = convert(Vector{Rational}, loaded_mps.c)
  objective = convert(Vector{fmpq}, objective)
  k = convert(Rational, loaded_mps.c0)
  k = convert(fmpq, k)
  int_var = Int64[]

  for (i, v_type) in enumerate(loaded_mps.vartypes)
      if QPSReader.VTYPE_Integer == v_type || QPSReader.VTYPE_Binary == v_type
          push!(int_var, i)
      end
  end

  return MixedIntegerLinearProgram{fmpq}(
    h_description,
    objective;
    integer_variables=int_var
  )
end

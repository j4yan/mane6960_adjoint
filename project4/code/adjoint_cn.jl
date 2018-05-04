export unsteadyAdjointSolveCN
export cmptTotalDerivOfFunc

function unsteadyAdjointSolveCN{T}(solver::EulerSolver{T},
                                   time::T, 
                                   nsteps::Int;
                                   savefile::AbstractString="solsave.dat")
    dt = time/nsteps

    fsave = open(savefile, "r")
    q_all = zeros(3, solver.sbp.numnodes, solver.numelem, nsteps+1)

    for i = 0 : nsteps
        # q = view(q_all, :, :, :, i + 1)
        q = zeros(Float64, 3, solver.sbp.numnodes, solver.numelem)
        read!(fsave, q)
        q_all[:,:,:, i+1] = q[:,:,:]
    end

    psi_all = zeros(Float64, size(q_all))
    psi_old = zeros(3, solver.sbp.numnodes, solver.numelem)
    wi = 1.

    for i = nsteps : -1 : 0
        if i == nsteps || i == 0
            wi = 0.5
        end
        psi = view(psi_all, :, :, :, i+1)
        q = view(q_all, :, :, :, i+1)
        stepCNAdj!(solver, q, psi, psi_old, wi, dt)
        psi_old = view(psi_all, :, :, :, i+1)
    end

    return psi_all
    return nothing  
end

function stepCNAdj!{T, Tsol, Tadj}(solver::EulerSolver{T}, 
                                   q::AbstractArray{Tsol, 3},
                                   psi::AbstractArray{Tadj, 3},
                                   psiold::AbstractArray{Tadj, 3},
                                   wn::Float64,
                                   dt::T)
    src = 0.0
    Jac = calcStateJacobian(solver, src, q)
    scale!(Jac, 0.5)
    Jac_rhs = deepcopy(Jac)
    addScaledMassMatrix!(solver, 1.0/dt, Jac)
    addScaledMassMatrix!(solver, -1.0/dt, Jac_rhs)
    rhs = zeros(Float64, size(psi))
    psiold_vec = reshape(psiold, length(psiold))
    rhs_vec = reshape(rhs, length(rhs))

    rhs[:]= calcdfuncdq(solver, q)
    rhs[:] = rhs[:] * wn
    rhs[:] = rhs[:] - Jac_rhs' * psiold_vec
    psi[:] = Jac' \ rhs_vec
    return nothing
end

"""
Compute the integrand in functional integral
    f = ∫ g dt
    g = 0.5 ∫ κ * (p - ptarg)^2 dx
"""
function calcFunctionalIntegrand{T, Tsol}(solver::EulerSolver{T},
                                          q::AbstractArray{Tsol, 3};
                                          sensor_x::T=0.85, 
                                          sensor_sig::T=0.05,
                                          press_targ::T=0.2)
    integrand = zeros(Tsol, solver.sbp.numnodes)
    scaled_int = zeros(integrand)
    fun_step = 0.0
    for k = 1 : solver.numelem
        for i = 1 : solver.sbp.numnodes
            dx = (solver.x[1,i,k] - sensor_x)/sensor_sig
            kern = exp(-dx*dx)
            dpress = calcPressure(view(q,:,i,k)) - press_targ
            integrand[i] = 0.5*dpress*dpress*kern*solver.jac[i,k]
        end
        fill!(scaled_int, zero(Tsol))
        volumeIntegrateElement!(solver.sbp, integrand, scaled_int)
        fun_step += sum(scaled_int)
    end
    return fun_step
end

"""
Compute the partial derivative of integrand, g, of functional integral w.r.t. q
    f = ∫ g dt
    g = 0.5 ∫ κ * (p - ptarg)^2 dx
"""
function calcdfuncdq{T, Tsol}(solver::EulerSolver{T},
                              q::AbstractArray{Tsol, 3};
                              sensor_x::T=0.85, 
                              sensor_sig::T=0.05,
                              press_targ::T=0.2)
    dfuncdq = zeros(Float64, size(q))
    qc = zeros(Complex128, size(q)) 
    qc[:] = q[:]
    eps = 1e-16

    for i = 1 : length(qc)
        qc[i] += Complex128(0, eps) 
        f = calcFunctionalIntegrand(solver, qc, 
                                    sensor_x=sensor_x, 
                                    sensor_sig=sensor_sig, 
                                    press_targ=press_targ)
        qc[i] -= Complex128(0, eps) 
        dfuncdq[i] = imag(f) / eps
    end
    return dfuncdq
end

function cmptTotalDerivOfFunc{T, Tsrc, Tsol, Tadj}(solver::EulerSolver{T},
                                                   dt::T,
                                                   src_all::AbstractArray{Tsrc, 1},
                                                   q_all::AbstractArray{Tsol, 4},
                                                   psi_all::AbstractArray{Tadj, 4})
    nsteps = size(q_all, 4) - 1
    dJds = zeros(nsteps + 1)
    @assert(nsteps == size(psi_all, 4) - 1)
    @assert(nsteps == length(src_all) - 1)
    # k = 0
    # src = src_all[k+1]
    # qcur = view(q_all, :, :, :, k+1)
    # psinew = view(psi_all, :, :, :, k+2)
    # psicur = zeros(size(psinew))
    # dJds[k] = cmptTotalDerivOfFunc(solver, dt, src, qcur, psicur, psinew)

    k = nsteps
    src = src_all[k+1]
    qcur = view(q_all, :, :, :, k+1)
    psicur = view(psi_all, :, :, :, k+1)
    psinew = zeros(size(psicur))
    dJds[k] = cmptTotalDerivOfFunc(solver, dt, src, qcur, psicur, psinew)

    for k = 1 : nsteps - 1
        src = src_all[k]
        qcur = view(q_all, :, :, :, k)
        qnew = view(q_all, :, :, :, k+1)
        psicur = view(psi_all, :, :, :, k)
        psinew = view(psi_all, :, :, :, k+1)
        dJds[k] = cmptTotalDerivOfFunc(solver, dt, src, qcur, psicur, psinew)
    end
    return dJds
end

"""
Compute the total derivative of the functional w.r.t. source.
"""
function cmptTotalDerivOfFunc{T, Tsrc, Tsol, Tadj}(solver::EulerSolver{T},
                                                   dt::T,
                                                   src::Tsrc,
                                                   qcur::AbstractArray{Tsol, 3},
                                                   psicur::AbstractArray{Tadj, 3},
                                                   psinew::AbstractArray{Tadj, 3})
    dRds = calcdRds(solver, src, qcur)
    dRds[:] *= dt
    dRds_vec = reshape(dRds, length(dRds))
    psicur_vec = reshape(psicur, length(psicur))
    psinew_vec = reshape(psinew, length(psinew))
    dJds = -dot(psicur_vec + psinew_vec, dRds_vec)
    return dJds
end

"""
Using complex-step to compute the partial derivative of spatial residual
with respect to the source. Note this derivative is not scaled by time step.
"""
function calcdRds{T, Tsrc, Tsol}(solver::EulerSolver{T}, 
                                 src::Tsrc,
                                 q::AbstractArray{Tsol, 3})
    eps = 1e-16
    srcc = Complex128(src, 0.0)
    srcc += Complex128(0., eps)
    res = zeros(Complex128, size(q))
    calcWeakResidual!(solver, srcc, q, res)
    dRds = zeros(size(q))
    for i = 1 : length(res)
        dRds[i] = imag(res[i]) / eps
    end
    return dRds 
end

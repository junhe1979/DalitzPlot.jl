using DalitzPlot,StaticArrays
using Test

@testset "DalitzPlot.jl" begin
    
    function amp_fp2KKp(kf, ch, para)

        ebm1 = sqrt(para.pW^2 + ch.p_i[1]^2)
        tecm = sqrt((ebm1 + ch.p_i[2])^2 - para.pW^2)
        #flux factor for cross section
        fac = 1 / (4 * para.pW * ch.p_i[2] * (2 * pi)^5)
        #fermion, aevrage on initial particele 
        fac = fac * 2.0 * ch.p_i[2] * 2.0 * ch.p_f[2] / 2.0
        #质心系到实验室系
        GAM = (ebm1 + ch.p_i[2]) / tecm
        ETA = para.pW / tecm
        for i in 1:3
            EI = kf[4, i]
            PX = kf[1, i]
            kf[1, i] = GAM * PX + ETA * EI
            kf[4, i] = ETA * PX + GAM * EI
            kf[5, i] = ch.p_f[i] #注意GENEV产生的kf[5]是三动量的大小
        end
    
        l = para.l::Float64
        p_p, p_f, p_K = ch.p_i[2], ch.p_i[1], ch.p_f[1]
        p_rho, p_omega, p_sigma = 0.77526, 0.78626, 0.5
        EB = 2.0 * p_K - p_f
        # k.k
        pf = SVector{5,Float64}([para.pW 0.0 0.0 ebm1 ch.p_i[1]]) #f0
        pp = SVector{5,Float64}([0.0 0.0 0.0 ch.p_i[2] ch.p_i[2]]) #p
        kKm = SVector{5,Float64}(kf[:, 1]) #K-
        kKp = SVector{5,Float64}(kf[:, 2]) #K+
        kp = SVector{5,Float64}(kf[:, 3]) #p
    
        q = pp - kp
        q2 = cdot(q, q)
        p = pf - kKm
        p2 = cdot(p, p)
    
        # propagator with monopole form factor
        gfKK = 5.0 #sqrt(64.0 * pi * p_f * sqrt(p_K * EB))
        #gfKK = gfKK / (p2 - p_K^2+im*p_K*5.3169e-17)
        gfKK = gfKK / (p2 - p_K^2 + im * p_K * 5.3169e-6)
    
        #aa = 1.0 / sqrt(p_K * EB)
        #gfKK = sqrt(8.0 / p_f) * p_K * sqrt(8.0 * pi / aa) / (kKm[1]^2 + kKm[2]^2 + kKm[3]^2 + 1.0 / aa^2)
    
        grho = 3.02 * (-3.1) / p_rho^2 * (p_rho^2 - l^2) / (q2 - l^2)
        gomega = -3.02 * 3.0 * (-3.1) / p_omega^2 * (p_omega^2 - l^2) / (q2 - l^2)
        gsigma = -3.65 * 12.8 / (q2 - p_sigma^2) * (p_sigma^2 - l^2) / (q2 - l^2)
        gq = (gfKK * (grho + gomega)) * q
        Mq = GS(ComplexF64.(gq))
        M = Mq + gfKK * gsigma * III
        #M = III #Mq * Mq 
        total0 = 0e0
        for ilNf in -1:2:1
            Uf = Ubc(kp, ilNf)
            for ilNi in -1:2:1
                Ui = Uc(pp, ilNi)
                total0 += abs2((Uf*M*Ui)[1])
            end
        end
    
        total = total0 * fac
    
        return total
    end
    p_K, p_p, p_f, p_L = 0.493677, 0.938272081, 0.987, 1.115683
    ch = (p_i=[p_f, p_p], p_f=[p_K, p_K, p_p],
        Lp_i=["f_0(980)", "p"], Lp_f=["K^-", "K^+", "p"],
        amp=amp_fp2KKp) #初末态粒子质量

    pW = 3.0
    ebm =  sqrt(pW^2 + ch.p_i[1]^2) #
    tecm = sqrt((ebm + ch.p_i[2])^2 - pW^2)

    res=Xsection(tecm, ch,nevtot=Int64(1e6),para=(pW=pW, l=1.0))
    @show pW, res.cs0
    plotD(res,ch,axes=[1, 3])

end

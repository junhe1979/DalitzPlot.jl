    amp_fp2KKp(kf, ch, para)=1.
    p_K, p_p, p_f, p_L = 0.493677, 0.938272081, 0.987, 1.115683
    ch = (p_i=[p_f, p_p], p_f=[p_K, p_K, p_p],
        Lp_i=["f_0(980)", "p"], Lp_f=["K^-", "K^+", "p"],
        amp=amp_fp2KKp) #初末态粒子质量

    pW = 3.0
    ebm =  sqrt(pW^2 + ch.p_i[1]^2) #
    tecm = sqrt((ebm + ch.p_i[2])^2 - pW^2)

    res=Xsection(tecm, ch,nevtot=Int64(1e6),para=(pW=pW, l=1.0))
    @show pW, res.cs0
    plotDP.plotD(res,ch,axes=[1, 3])
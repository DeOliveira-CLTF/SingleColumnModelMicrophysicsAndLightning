! ============================================================================
! mpas_atmphys_lightning.F90
!
! Parametrizacao diagnostica de taxa de relampago (flash rate) para o
! MONAN/MPAS-A. Segue o padrao de modulo de fisica do MPAS (uma subrotina
! "driver_*" que decide o esquema por coluna, chamada apos a microfisica
! e apos o esquema de cumulus no driver de fisica).
!
! Referencias:
!   McCaul, E.W. et al. (2009), WAF, 24, 709-729.  (forma original: linear,
!     fluxo w*qg no nivel de -15C -- ver flashrate_mccaul_column, mantida
!     so para comparacao)
!   Price, C. & Rind, D. (1992), JGR, 97, 9919-9933. (via altura do topo)
!   Price, C. & Rind, D. (1993), JAM, 32, 170-181.  (particao IC/CG)
!   Barthe, C. et al. (2010), JGR, 115, D24202. (recalibracao Price-Rind
!     via regressao contra altura do topo, r=0.74, coef_alpha_pr/coef_beta_pr)
!   Thompson, G. et al. -- esquema de microfisica (N0,g diagnosticado a
!     partir de qg, usado no fechamento de tamanho do graupel no EP)

!
! Paulo Kubota / INPE-CPTEC - esqueleto para uso didatico (MET-576-4) e
! prototipo de integracao operacional no MONAN.
! ============================================================================
module mpas_atmphys_lightning

   implicit none
   private

   public :: driver_lightning
   public :: flashrate_mccaul_column
   public :: flashrate_pricerind
   public :: calc_ic_cg_ratio
   public :: flashrate_ep_column

   integer, parameter :: rk = selected_real_kind(12)

   ! -- coeficientes empiricos (expor como namelist na integracao real) -----
   real(rk), parameter :: coef_a_mccaul   = 0.042_rk    ! McCaul et al. (2009), eq. F1: slope
                                                         ! da regressao linear FR = k1*(w*qg)|_-15C

   real(rk), parameter :: coef_alpha_pr   = 1.71e-9_rk  ! FR_total = alpha * H^beta (H em km)
   real(rk), parameter :: coef_beta_pr    = 8.7_rk
   ! Recalibrado com Barthe et al. (2010, JGR) -- regressao entre altura do
   ! topo (cloud top height) e a TAXA DE FLASH TOTAL observada/simulada
   ! (r=0.74), a partir de simulacoes WRF de 2 casos (STERAO-A e Alabama do
   ! Norte). Substitui os coeficientes originais de Price & Rind (1992)
   ! (alpha=3.44e-5, beta=4.9). ATENCAO: calibracao extratropical (EUA),
   ! nao tropical -- ainda precisa de recalibracao regional com RINDAT/GLM
   ! sobre casos de MCS na Bacia do Prata / Amazonia antes de uso
   ! operacional no MONAN, como ja apontado para os demais coeficientes.

   real(rk), parameter :: t_peak_nic      = -15.0_rk    ! nivel-alvo do McCaul original (2009):
                                                         ! pico de eficiencia NIC graupel-cristal
                                                         ! (Takahashi 1978/Saunders), onde o paper
                                                         ! avalia o fluxo w*qg (nao integrado)

   ! -- transicao scale-aware (sigmoide em dx) ------------------------------
   ! Em vez de chaveamento rigido if(dx_km > dx_threshold_km), a fracao
   ! "convectiva parametrizada" sigma_cu cresce suavemente de 0 (malha fina,
   ! convection-permitting) a 1 (malha grossa, cumulus dominante). Evita
   ! salto artificial na faixa de resolucao "cinzenta" (~3-10 km) onde o
   ! Grell-Freitas ja faz sua propria transicao suave (fallback ate que o
   ! sigma real do GF seja exposto no pool de diagnostico).
   real(rk), parameter :: dx_threshold_km  = 6.0_rk     ! centro da transicao (km)
   real(rk), parameter :: dx_transition_km = 2.0_rk     ! largura da transicao (km)

   ! -- Potencial de Eletrificacao (EP) -------------------------------------
   ! Proxy fisico que substitui o McCaul (massa simples) na branch de grade
   ! explicita/convection-permitting, portado do prototipo Python
   ! (ep_prototype.py) ja validado e com os 2 bugs corrigidos:
   !   (1) conversao de EW (agua liquida efetiva) sem fator 1e-3 espurio
   !   (2) magnitude da carga transferida cai perto da temperatura de
   !       reversao (dist_factor), e cresce com a intensidade do updraft
   ! Teoria: separacao de carga nao-indutiva (NIC) graupel-cristal
   ! (Takahashi 1978; Saunders & Peck 1998; Jayaratne et al. 1983).
   !
   ! ATUALIZACAO (curso MET-756-4, apostila do esquema Thompson): a
   ! concentracao numerica e o diametro caracteristico do graupel NAO usam
   ! mais um diametro fixo assumido -- sao derivados do parametro de
   ! interceptacao N0,g diagnosticado a partir do proprio qg, exatamente
   ! como o Thompson scheme faz internamente:
   !   N0,g = max(1e4, min(200/qg, 5e6))  [m^-4], qg em kg/kg
   ! seguido do fechamento da distribuicao exponencial (rho_graupel=400
   ! kg/m3, igual ao valor usado no Thompson) para obter lambda_g, Ng e o
   ! diametro caracteristico Dg=1/lambda_g. Isso deixa Ng/Vg consistentes
   ! com a microfisica que o MONAN de fato roda, em vez de uma aproximacao
   ! paralela com diametro fixo.
   !
   ! ATENCAO: constantes de microfisica assumida (D0_ICE, densidades) e
   ! coeficientes da lei de potencia (coef_a_ep, coef_b_ep) ainda sao
   ! ORDENS DE GRANDEZA DA LITERATURA, NAO CALIBRADAS -- precisam de ajuste
   ! contra dados reais (RINDAT/GLM) antes de uso operacional.
   real(rk), parameter :: rho_graupel_ep  = 400.0_rk    ! kg/m3, densidade do graupel (Thompson)
   real(rk), parameter :: d0_ice_ep       = 2.0e-4_rk   ! m, diametro medio assumido cristal
   real(rk), parameter :: rho_ice_ep      = 900.0_rk    ! kg/m3, densidade assumida do gelo
   real(rk), parameter :: v_ice_ep        = 0.5_rk      ! m/s, veloc. terminal tipica cristal
   real(rk), parameter :: d0_graupel_fallback_ep = 2.0e-3_rk  ! m, usado so quando qg<=0

   real(rk), parameter :: t_top_ep        = -40.0_rk    ! limite superior camada fase mista (C)
   real(rk), parameter :: t_bot_ep        = 0.0_rk      ! limite inferior camada fase mista (C)

   ! lei de potencia FR_EP = a_ep * EP^b_ep -- PLACEHOLDER, calibrar depois
   real(rk), parameter :: coef_a_ep       = 1.0e-3_rk
   real(rk), parameter :: coef_b_ep       = 0.5_rk

contains

   ! ------------------------------------------------------------------------
   ! driver_lightning: decide o esquema por coluna com base na resolucao
   ! local (dx) e no esquema de cumulus ativo (cu_active).
   !
   ! Se cu_active = .true. (convecção parametrizada, sem graupel de grade
   ! confiavel) -> usa Price-Rind (precisa da altura do topo convectivo,
   ! ja diagnosticada pelo esquema de cumulus, ex. Grell-Freitas: cldtop_km).
   !
   ! Se cu_active = .false. (convection-permitting, graupel resolvido)
   ! -> usa McCaul via integral de massa de gelo na regiao de carga.
   ! ------------------------------------------------------------------------
   subroutine driver_lightning(nz, z_m, t_c, rho, w, qg, qs, qc, dx_km, &
                                cu_active, cldtop_km, &
                                flash_rate, flash_rate_cg, flash_rate_ic)
      integer,  intent(in)  :: nz
      real(rk), intent(in)  :: z_m(nz)      ! altura geometrica (m), niveis do modelo
      real(rk), intent(in)  :: t_c(nz)      ! temperatura (C)
      real(rk), intent(in)  :: rho(nz)      ! densidade do ar (kg/m3)
      real(rk), intent(in)  :: w(nz)        ! velocidade vertical (m/s)
      real(rk), intent(in)  :: qg(nz)       ! razao de mistura graupel (kg/kg)
      real(rk), intent(in)  :: qs(nz)       ! razao de mistura neve (kg/kg)
      real(rk), intent(in)  :: qc(nz)       ! razao de mistura agua liquida de nuvem (kg/kg)
      real(rk), intent(in)  :: dx_km        ! resolucao local da malha (km)
      logical,  intent(in)  :: cu_active    ! esquema de cumulus ativo nesta coluna?
      real(rk), intent(in)  :: cldtop_km    ! altura do topo convectivo (km), se cu_active

      real(rk), intent(out) :: flash_rate     ! flashes / min (total, IC+CG)
      real(rk), intent(out) :: flash_rate_cg  ! flashes / min (nuvem-solo)
      real(rk), intent(out) :: flash_rate_ic  ! flashes / min (intra-nuvem)

      real(rk) :: ratio_ic_cg, h_eff, sigma_cu
      real(rk) :: fr_pr, fr_pr_cg, fr_pr_ic, fr_pr_total
      real(rk) :: fr_ep, fr_ep_cg, fr_ep_ic
      real(rk) :: ep_mag, ep_signed

      ! fracao "convectiva parametrizada" (0 = malha fina/convection-
      ! permitting -> so EP; 1 = malha grossa/cumulus dominante ->
      ! so Price-Rind). Se cu_active = .false. em toda a coluna (esquema
      ! de cumulus desligado no dominio), forca sigma_cu=0 independente de dx.
      if (.not. cu_active) then
         sigma_cu = 0.0_rk
      else
         sigma_cu = sigmoid_dx(dx_km)
      end if

      ! --- via Price-Rind/Barthe (altura do topo convectivo) -------------
      ! CORRECAO: flashrate_pricerind retorna o flash rate TOTAL (assim
      ! como a regressao de Barthe et al. 2010 contra cloud top height, e
      ! a formula original de Price & Rind 1992) -- nao apenas a fracao CG
      ! como uma versao anterior deste modulo assumia por engano. O total
      ! e dividido entre CG e IC usando a razao IC/CG (Price & Rind 1993).
      h_eff = cldtop_km
      call flashrate_pricerind(h_eff, fr_pr_total)
      call calc_ic_cg_ratio(h_eff, ratio_ic_cg)
      fr_pr_cg = fr_pr_total / (1.0_rk + ratio_ic_cg)
      fr_pr_ic = fr_pr_total - fr_pr_cg
      fr_pr    = fr_pr_total

      ! --- via EP (potencial de eletrificacao, grade explicita) -----------
      ! Substitui o proxy de massa simples do McCaul nesta branch: o EP
      ! agrega sensibilidade a fase de crescimento do graupel (seco/umido,
      ! via reversao de sinal dependente de LWC/qc) que o McCaul nao ve.
      ! flashrate_mccaul_column permanece disponivel (publica) para
      ! comparacao/diagnostico, mas nao e mais usada pelo driver.
      call flashrate_ep_column(nz, z_m, t_c, rho, w, qg, qs, qc, ep_mag, ep_signed, fr_ep)
      ! sem info explicita de topo, usa particao IC/CG generica (~90/10)
      ! ate haver calibracao regional com RINDAT/GLM
      fr_ep_ic = fr_ep * 0.90_rk
      fr_ep_cg = fr_ep * 0.10_rk

      ! --- mistura ponderada pela fracao convectiva (scale-aware) ---------
      flash_rate    = sigma_cu * fr_pr    + (1.0_rk - sigma_cu) * fr_ep
      flash_rate_cg = sigma_cu * fr_pr_cg + (1.0_rk - sigma_cu) * fr_ep_cg
      flash_rate_ic = sigma_cu * fr_pr_ic + (1.0_rk - sigma_cu) * fr_ep_ic

   end subroutine driver_lightning

   ! ------------------------------------------------------------------------
   ! sigmoid_dx: transicao suave 0->1 centrada em dx_threshold_km, com
   ! largura dx_transition_km. Fallback ate que a fracao de area convectiva
   ! scale-aware real do Grell-Freitas esteja disponivel no pool de
   ! diagnostico do MONAN.
   ! ------------------------------------------------------------------------
   pure function sigmoid_dx(dx_km) result(sigma)
      real(rk), intent(in) :: dx_km
      real(rk) :: sigma
      sigma = 1.0_rk / (1.0_rk + exp(-(dx_km - dx_threshold_km) / dx_transition_km))
   end function sigmoid_dx

   ! ------------------------------------------------------------------------
   ! flashrate_ep_column: Potencial de Eletrificacao (EP), portado do
   ! prototipo Python ep_prototype.py (ja validado nos 3 experimentos:
   ! reversao de sinal com LWC, curva Trev(EW), sensibilidade monotonica a
   ! w/qg). Substitui o McCaul como proxy fisico da branch de grade
   ! explicita (convection-permitting).
   !
   ! Teoria (NIC graupel-cristal, Takahashi 1978 / Saunders & Peck 1998):
   ! a taxa de separacao de carga depende da taxa de colisao
   ! (Ng*Ni*|Vg-Vi|*secao_choque), da eficiencia de colisao Eci(T) (pico em
   ! ~-15C), e do SINAL/magnitude da carga transferida, que depende de quao
   ! longe T esta da temperatura de reversao Trev(EW) -- que por sua vez
   ! desloca com o conteudo de agua liquida efetiva (aqui aproximado por qc).
   !
   ! Unidades de entrada: qg, qs, qc em kg/kg (convertidas internamente para
   ! g/kg, como no prototipo Python). Saida: EP (magnitude, unidade
   ! arbitraria ate calibracao) e EP_signed (guardado para uso futuro em
   ! polaridade CG+/CG- e LNOx), alem do flash_rate = a_ep * EP^b_ep.
   ! ------------------------------------------------------------------------
   subroutine flashrate_ep_column(nz, z_m, t_c, rho, w, qg, qs, qc, &
                                   ep_mag, ep_signed, flash_rate)
      integer,  intent(in)  :: nz
      real(rk), intent(in)  :: z_m(nz), t_c(nz), rho(nz), w(nz)
      real(rk), intent(in)  :: qg(nz), qs(nz), qc(nz)
      real(rk), intent(out) :: ep_mag       ! potencial de eletrificacao (magnitude)
      real(rk), intent(out) :: ep_signed    ! idem, com polaridade (uso futuro CG+/-)
      real(rk), intent(out) :: flash_rate

      real(rk) :: m_ice, xsec_k, xsec_km1
      real(rk) :: qs_gkg, qc_gkg
      real(rk) :: ng, ni, vg, dv, ew, trev, eci, dist_factor, updraft_factor
      real(rk) :: dg_char_k, dg_char_km1
      real(rk) :: dqdt_k, dqdt_km1, sign_k, sign_km1, dz
      logical  :: in_layer_k, in_layer_km1
      integer  :: k

      m_ice     = (3.14159265358979_rk/6.0_rk) * rho_ice_ep     * d0_ice_ep**3

      ep_mag    = 0.0_rk
      ep_signed = 0.0_rk

      do k = 2, nz
         in_layer_k   = (t_c(k)   <= t_bot_ep) .and. (t_c(k)   >= t_top_ep)
         in_layer_km1 = (t_c(k-1) <= t_bot_ep) .and. (t_c(k-1) >= t_top_ep)

         if (in_layer_k .or. in_layer_km1) then

            ! --- nivel k -----------------------------------------------
            call graupel_size_ep(qg(k), rho(k), ng, dg_char_k)
            qs_gkg = qs(k) * 1000.0_rk
            qc_gkg = qc(k) * 1000.0_rk

            ni  = rho(k) * (qs_gkg * 1.0e-3_rk) / m_ice
            vg  = fallspeed_graupel_ep(dg_char_k * 1000.0_rk)
            dv  = abs(vg - v_ice_ep)
            xsec_k = 3.14159265358979_rk * (0.5_rk*dg_char_k + 0.5_rk*d0_ice_ep)**2
            ew  = max(qc_gkg, 0.0_rk) * rho(k)
            trev = reversal_temperature_ep(ew)
            eci  = collision_efficiency_ep(t_c(k))
            dist_factor   = tanh(abs(t_c(k) - trev) / 5.0_rk)
            updraft_factor = 1.0_rk + max(w(k), 0.0_rk) / 5.0_rk

            dqdt_k = ng * ni * dv * xsec_k * eci * dist_factor * updraft_factor
            if (t_c(k) > trev) then
               sign_k = -1.0_rk
            else
               sign_k = 1.0_rk
            end if

            ! --- nivel k-1 -----------------------------------------------
            call graupel_size_ep(qg(k-1), rho(k-1), ng, dg_char_km1)
            qs_gkg = qs(k-1) * 1000.0_rk
            qc_gkg = qc(k-1) * 1000.0_rk

            ni  = rho(k-1) * (qs_gkg * 1.0e-3_rk) / m_ice
            vg  = fallspeed_graupel_ep(dg_char_km1 * 1000.0_rk)
            dv  = abs(vg - v_ice_ep)
            xsec_km1 = 3.14159265358979_rk * (0.5_rk*dg_char_km1 + 0.5_rk*d0_ice_ep)**2
            ew  = max(qc_gkg, 0.0_rk) * rho(k-1)
            trev = reversal_temperature_ep(ew)
            eci  = collision_efficiency_ep(t_c(k-1))
            dist_factor   = tanh(abs(t_c(k-1) - trev) / 5.0_rk)
            updraft_factor = 1.0_rk + max(w(k-1), 0.0_rk) / 5.0_rk

            dqdt_km1 = ng * ni * dv * xsec_km1 * eci * dist_factor * updraft_factor
            if (t_c(k-1) > trev) then
               sign_km1 = -1.0_rk
            else
               sign_km1 = 1.0_rk
            end if

            dz = z_m(k) - z_m(k-1)

            ep_mag    = ep_mag    + 0.5_rk*(abs(dqdt_k) + abs(dqdt_km1))*dz
            ep_signed = ep_signed + 0.5_rk*(dqdt_k*sign_k + dqdt_km1*sign_km1)*dz
         end if
      end do

      if (ep_mag > 0.0_rk) then
         flash_rate = coef_a_ep * ep_mag**coef_b_ep
      else
         flash_rate = 0.0_rk
      end if

   end subroutine flashrate_ep_column

   ! ------------------------------------------------------------------------
   ! n0g_thompson: parametro de interceptacao do graupel diagnosticado a
   ! partir de qg (kg/kg), exatamente como no Thompson scheme (curso
   ! MET-756-4, apostila "Gregory Thompson Scheme"):
   !   N0,g = max(1e4, min(200/qg, 5e6))  [m^-4]
   ! Correntes ascendentes fortes (mais qg) -> N0,g menor -> espectro desloca
   ! para particulas maiores -> queda mais rapida (sem precisar aumentar N0,g
   ! como um parametro fixo faria).
   ! ------------------------------------------------------------------------
   pure function n0g_thompson(qg_kgkg) result(n0g)
      real(rk), intent(in) :: qg_kgkg
      real(rk) :: n0g, qg_safe
      qg_safe = max(qg_kgkg, 1.0e-12_rk)
      n0g = max(1.0e4_rk, min(200.0_rk/qg_safe, 5.0e6_rk))
   end function n0g_thompson

   ! ------------------------------------------------------------------------
   ! graupel_size_ep: concentracao numerica (Ng) e diametro caracteristico
   ! (Dg=1/lambda_g) do graupel, a partir do fechamento da distribuicao
   ! exponencial N(D)=N0,g*exp(-lambda_g*D) com N0,g diagnosticado por
   ! n0g_thompson(qg) e rho_graupel_ep=400 kg/m3 (igual ao Thompson).
   ! lambda_g^4 = 6*pi*rho_graupel*N0,g / (rho_ar * qg)
   ! ------------------------------------------------------------------------
   subroutine graupel_size_ep(qg_kgkg, rho_k, ng, dg_char)
      real(rk), intent(in)  :: qg_kgkg, rho_k
      real(rk), intent(out) :: ng, dg_char
      real(rk) :: n0g, lambda_g

      if (qg_kgkg <= 1.0e-12_rk) then
         ng = 0.0_rk
         dg_char = d0_graupel_fallback_ep  ! sem massa de graupel -> Ng=0 zera a colisao de qualquer forma
         return
      end if

      n0g = n0g_thompson(qg_kgkg)
      lambda_g = (6.0_rk*3.14159265358979_rk*rho_graupel_ep*n0g / (rho_k*qg_kgkg))**0.25_rk
      ng = n0g / lambda_g
      dg_char = 1.0_rk / lambda_g
   end subroutine graupel_size_ep

   ! ------------------------------------------------------------------------
   ! fallspeed_graupel_ep: velocidade terminal do graupel (m/s), crescendo
   ! com o diametro caracteristico Dg (mm) -- agora derivado fisicamente de
   ! N0,g(qg) (graupel_size_ep) em vez de um diametro fixo assumido. Forma
   ! de lei de potencia, calibrar depois.
   ! ------------------------------------------------------------------------
   pure function fallspeed_graupel_ep(dg_mm) result(vg)
      real(rk), intent(in) :: dg_mm
      real(rk) :: vg
      vg = 3.0_rk * max(dg_mm, 1.0e-3_rk)**0.15_rk + 1.5_rk
   end function fallspeed_graupel_ep

   ! ------------------------------------------------------------------------
   ! collision_efficiency_ep: eficiencia de colisao/riming graupel-cristal,
   ! pico proximo a -15C (Takahashi 1978), decaindo para os dois lados.
   ! ------------------------------------------------------------------------
   pure function collision_efficiency_ep(t_c_k) result(eci)
      real(rk), intent(in) :: t_c_k
      real(rk) :: eci
      eci = exp(-0.5_rk * ((t_c_k + 15.0_rk) / 8.0_rk)**2)
   end function collision_efficiency_ep

   ! ------------------------------------------------------------------------
   ! reversal_temperature_ep: temperatura de reversao de sinal de carga (C)
   ! em funcao do conteudo de agua liquida efetivo (EW, g/m3), forma
   ! simplificada inspirada em Saunders & Peck (1998). Clampada em
   ! [-20, -3] C (faixa fisica tipica).
   ! ------------------------------------------------------------------------
   pure function reversal_temperature_ep(ew_gm3) result(trev)
      real(rk), intent(in) :: ew_gm3
      real(rk) :: trev
      trev = -15.0_rk + 6.0_rk * ew_gm3
      trev = min(max(trev, -20.0_rk), -3.0_rk)
   end function reversal_temperature_ep

   ! ------------------------------------------------------------------------
   ! flashrate_mccaul_column: forma LITERAL de McCaul et al. (2009), eq. F1
   ! -- descoberta ao ler o paper original (curso MET-756-4): a relacao
   ! publicada e LINEAR (FR = k1*(w*qg)), nao lei de potencia, e o fluxo
   ! de graupel e avaliado num UNICO NIVEL, exatamente a -15C (nao
   ! integrado em -5/-40C como uma versao anterior deste modulo fazia --
   ! aquela era uma extensao fisica minha, nao o McCaul original).
   ! Mantida publica so para comparacao/diagnostico; o driver usa o EP.
   ! ------------------------------------------------------------------------
   subroutine flashrate_mccaul_column(nz, z_m, t_c, rho, w, qg, qs, flash_rate)
      integer,  intent(in)  :: nz
      real(rk), intent(in)  :: z_m(nz), t_c(nz), rho(nz), w(nz), qg(nz), qs(nz)
      real(rk), intent(out) :: flash_rate

      real(rk) :: w_m15, qgs_m15, frac, graupel_flux_gkg
      integer  :: k, k_lo
      logical  :: found

      ! localiza o par de niveis que envolve T=-15C (busca a mudanca de
      ! sinal de (T - (-15)) entre niveis consecutivos)
      found = .false.
      k_lo = 1
      do k = 2, nz
         if ((t_c(k) - t_peak_nic) * (t_c(k-1) - t_peak_nic) <= 0.0_rk) then
            k_lo = k - 1
            found = .true.
            exit
         end if
      end do

      if (found) then
         if (abs(t_c(k) - t_c(k_lo)) > 1.0e-6_rk) then
            frac = (t_peak_nic - t_c(k_lo)) / (t_c(k) - t_c(k_lo))
         else
            frac = 0.5_rk
         end if
         frac = min(max(frac, 0.0_rk), 1.0_rk)

         w_m15   = w(k_lo)            + frac*(w(k)            - w(k_lo))
         qgs_m15 = (qg(k_lo)+qs(k_lo)) + frac*((qg(k)+qs(k)) - (qg(k_lo)+qs(k_lo)))
      else
         ! coluna nao atravessa -15C (rasa/quente ou muito fria em toda a
         ! extensao) -- usa o nivel mais proximo do alvo como aproximacao
         do k = 1, nz
            if (abs(t_c(k) - t_peak_nic) < abs(t_c(k_lo) - t_peak_nic)) k_lo = k
         end do
         w_m15   = w(k_lo)
         qgs_m15 = qg(k_lo) + qs(k_lo)
      end if

      ! graupel flux = w * qg no nivel de -15C, qg convertido para g/kg
      ! (convencao do paper original -- fluxos reportados na ordem de
      ! centenas, ex. "graupel flux of 400 m s-1", consistente com qg em
      ! g/kg e nao em kg/kg)
      graupel_flux_gkg = max(w_m15, 0.0_rk) * qgs_m15 * 1000.0_rk

      flash_rate = coef_a_mccaul * graupel_flux_gkg

   end subroutine flashrate_mccaul_column

   ! ------------------------------------------------------------------------
   ! flashrate_pricerind: FR_total = alpha * H^beta, H em km (topo
   ! convectivo). Coeficientes recalibrados com Barthe et al. (2010) --
   ! retorna o flash rate TOTAL (IC+CG), nao apenas CG (ver nota em
   ! driver_lightning sobre a correcao desta semantica).
   ! ------------------------------------------------------------------------
   subroutine flashrate_pricerind(h_km, flash_rate_total)
      real(rk), intent(in)  :: h_km
      real(rk), intent(out) :: flash_rate_total

      if (h_km > 0.0_rk) then
         flash_rate_total = coef_alpha_pr * h_km**coef_beta_pr
      else
         flash_rate_total = 0.0_rk
      end if
   end subroutine flashrate_pricerind

   ! ------------------------------------------------------------------------
   ! calc_ic_cg_ratio: razao IC/CG em funcao da altura do topo (Price & Rind
   ! 1993), valido aproximadamente para H entre 5.5 e 14 km; fora disso,
   ! clampado.
   ! ------------------------------------------------------------------------
   subroutine calc_ic_cg_ratio(h_km, ratio)
      real(rk), intent(in)  :: h_km
      real(rk), intent(out) :: ratio
      real(rk) :: h

      h = min(max(h_km, 5.5_rk), 14.0_rk)

      ratio = 0.021_rk*h**4 - 0.648_rk*h**3 + 7.49_rk*h**2 &
              - 36.54_rk*h + 63.09_rk

      ratio = max(ratio, 0.1_rk)  ! evita razao negativa/instavel nos extremos
   end subroutine calc_ic_cg_ratio

end module mpas_atmphys_lightning

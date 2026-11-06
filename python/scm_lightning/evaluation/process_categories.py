GROUPS={
 'Condensação/evap. de nuvem':['th_prw_vcd','ws_pcond'],
 'Deposição/sublimação':['th_pri_ide','th_prs_ide','th_prs_sde','th_prg_gde','ws_pidep','ws_psdep','ws_pgdep'],
 'Nucleação de gelo':['th_pri_inu','th_pri_iha','ws_pigen'],
 'Congelamento':['th_pri_wfz','th_pri_rfz','th_prg_rfz','th_homogeneous_freezing','ws_pihmf','ws_pfrzdtc','ws_pfrzdtr'],
 'Riming (água de nuvem)':['th_prs_scw','th_prg_scw','th_prg_gcw','ws_paacw_snow','ws_paacw_graupel'],
 'Coleta chuva–gelo':['th_collection_rci','th_collection_rcs','th_collection_rcg','ws_piacr','ws_praci','ws_psacr','ws_pracs','ws_pgacr'],
 'Coleta gelo–gelo':['th_prs_sci','ws_psaci','ws_pgaci','ws_pgacs'],
 'Autoconversão/agregação':['th_prr_wau','th_prs_iau','ws_praut','ws_psaut','ws_pgaut'],
 'Coleta quente':['th_prr_rcw','ws_pracw','ws_paacw_warm'],
 'Fusão':['th_prr_sml','th_prr_gml','th_instant_melting','ws_pimlt','ws_psmlt','ws_pgmlt','ws_pseml','ws_pgeml'],
 'Evaporação precipitação':['th_prv_rev','ws_prevp','ws_psevp','ws_pgevp'],
 'Gelo secundário':['th_prs_ihm','th_prg_ihm'],
 'Sedimentação':['sedimentation'],
 'Ajustes/limiares':['tendency_adjustments','positivity_clipping','final_clipping'],
 'Resíduo não atribuído':['unattributed_residual'],
}

CATEGORY={p:g for g,ps in GROUPS.items() for p in ps}

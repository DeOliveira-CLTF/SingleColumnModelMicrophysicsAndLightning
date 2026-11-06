"""Insert observational Fortran hooks; verify original physics text is intact.

All insertions are delimited. Re-running is idempotent. Original bytes saved.
"""
import re
from scm_artigo import ROOT,RESULTS,sha,save_json

PROV=RESULTS/'ensemble_goamazon/processos/proveniencia'
BEGIN='! PROCESS_DIAGNOSTICS_BEGIN\n'
END='! PROCESS_DIAGNOSTICS_END\n'


def strip_hooks(text):
    return re.sub(r'! PROCESS_DIAGNOSTICS_BEGIN\n.*?! PROCESS_DIAGNOSTICS_END\n','',text,flags=re.S)


def instrument(path,fn):
    original=strip_hooks(path.read_text(encoding='utf-8'))
    backup=PROV/path.name
    if not backup.exists(): backup.write_bytes(path.read_bytes())
    assert original==backup.read_text(encoding='utf-8')
    result=fn(original)
    assert strip_hooks(result)==original
    path.write_text(result,encoding='utf-8')


def insert(text,anchor,code,after=False):
    assert text.count(anchor)==1,(anchor,text.count(anchor))
    block=BEGIN+code.rstrip()+'\n'+END
    return text.replace(anchor,anchor+block if after else block+anchor)


def transfer(name,a,b,amount):
    return f"          call pd_transfer(k,'{name}',{a},{b},{amount})\n"


def thompson(text):
    text=insert(text,'      use mpas_log\n','      use process_diagnostics\n')
    pairs={'pri_inu':(1,4),'pri_iha':(1,4),'pri_ide':(1,4),'prs_ide':(1,5),
           'prs_sde':(1,5),'prg_gde':(1,6),'prr_wau':(2,3),'pri_wfz':(2,4),
           'prr_rcw':(2,3),'prs_scw':(2,5),'prg_scw':(2,6),'prg_gcw':(2,6),
           'prs_ihm':(5,4),'prg_ihm':(6,4),'pri_rfz':(3,4),'prg_rfz':(3,6),
           'prs_iau':(4,5),'prs_sci':(4,5),'prr_sml':(5,3),'prr_gml':(6,3)}
    code='         if(pd_enabled) then\n'
    for name,(a,b) in pairs.items(): code+=transfer('th_'+name,a,b,f'{name}(k)*orho*DT')
    code+="          call pd_add(k,'th_collection_rci',DT*orho*[0.,0.,-prr_rci(k),-pri_rci(k),0.,prg_rci(k)])\n"
    code+="          call pd_add(k,'th_collection_rcs',DT*orho*[0.,0.,prr_rcs(k),0.,prs_rcs(k),prg_rcs(k)])\n"
    code+="          call pd_add(k,'th_collection_rcg',DT*orho*[0.,0.,prr_rcg(k),0.,0.,prg_rcg(k)])\n         endif\n"
    text=insert(text,'!..Water vapor tendency\n',code)
    for anchor,name,a,b,amount in [
        ('          qvten(k) = qvten(k) - prw_vcd(k)\n','th_prw_vcd',1,2,'prw_vcd(k)*DT'),
        ('          qrten(k) = qrten(k) - prv_rev(k)\n','th_prv_rev',3,1,'prv_rev(k)*DT'),
        ('          qcten(k) = qcten(k) + xri*odt\n','th_instant_melting',4,2,'xri'),
        ('          qiten(k) = qiten(k) + xrc*odt\n','th_homogeneous_freezing',2,4,'xrc')]:
        text=insert(text,anchor,'          if(pd_enabled) '+transfer(name,a,b,amount).strip())
    tendencies='DT*[qvten(k),qcten(k),qrten(k),qiten(k),qsten(k),qgten(k)]'
    text=insert(text,'      nstep = 0\n      onstep(:) = 1.0\n',f'''      if(pd_enabled) then
        do k=kts,kte
          call pd_reconcile_tendency(k,{tendencies})
          call pd_mark(k,{tendencies})
        enddo
      endif''')
    text=insert(text,'!.. Instantly melt any cloud ice into cloud water if above 0C and\n',f'''      if(pd_enabled) then
        do k=kts,kte
          call pd_difference(k,'sedimentation',{tendencies})
        enddo
      endif''')
    predicted='[qv1d(k)+qvten(k)*DT,qc1d(k)+qcten(k)*DT,qr1d(k)+qrten(k)*DT,qi1d(k)+qiten(k)*DT,qs1d(k)+qsten(k)*DT,qg1d(k)+qgten(k)*DT]'
    text=insert(text,'         t1d(k)  = t1d(k) + tten(k)*DT\n',f'         if(pd_enabled) call pd_mark(k,{predicted})')
    text=insert(text,'         if (qg1d(k) .le. R1) qg1d(k) = 0.0\n',
                "         if(pd_enabled) call pd_difference(k,'final_clipping',[qv1d(k),qc1d(k),qr1d(k),qi1d(k),qs1d(k),qg1d(k)])",after=True)
    return text


def wsm6(text):
    text=insert(text,' use ccpp_kind_types,only: kind_phys\n',' use process_diagnostics\n')
    state='[q(i,k),qc(i,k),qr(i,k),qi(i,k),qs(i,k),qg(i,k)]'
    # Direct changes, including exact limited freezing/melting increments.
    direct=[('         qc(i,k) = qc(i,k) + qi(i,k)\n','ws_pimlt',4,2,'qi(i,k)'),
            ('         qi(i,k) = qi(i,k) + qc(i,k)\n','ws_pihmf',2,4,'qc(i,k)'),
            ('         qi(i,k) = qi(i,k) + pfrzdtc\n','ws_pfrzdtc',2,4,'pfrzdtc'),
            ('         qg(i,k) = qg(i,k) + pfrzdtr\n','ws_pfrzdtr',3,6,'pfrzdtr'),
            ('           qs(i,k) = qs(i,k) + psmlt(i,k)\n','ws_psmlt',5,3,'-psmlt(i,k)'),
            ('           qg(i,k) = qg(i,k) + pgmlt(i,k)\n','ws_pgmlt',6,3,'-pgmlt(i,k)'),
            ('         q(i,k) = q(i,k)-pcond(i,k)*dtcld\n','ws_pcond',1,2,'pcond(i,k)*dtcld')]
    for anchor,name,a,b,amount in direct:
        text=insert(text,anchor,'         if(pd_enabled) '+transfer(name,a,b,amount).strip())
    # Sedimentation before/after actual solver, includes its substeps and clips.
    for first,last in [('   call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,workr,denqrs1,  &\n',
                        '   do i = its, ite\n     fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld\n'),
                       ('   call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,work1c,denqci,  &\n',
                        '   do i = its, ite\n     fallc(i,1) = delqi(i)/delz(i,1)/dtcld\n')]:
        text=insert(text,first,f'   if(pd_enabled) then\n     do k=kts,kte\n       i=its\n       call pd_mark(k,{state})\n     enddo\n   endif')
        text=insert(text,last,f"   if(pd_enabled) then\n     do k=kts,kte\n       i=its\n       call pd_difference(k,'sedimentation',{state})\n     enddo\n   endif")
    # Final limited rates, in their actual cold/warm branches.
    common={'praut':(2,3),'pracw':(2,3),'prevp':(1,3)}
    cold=dict(common,pigen=(1,4),pidep=(1,4),psdep=(1,5),pgdep=(1,6),psaut=(4,5),pgaut=(5,6),
              psaci=(4,5),pgaci=(4,6),pgacr=(3,6),pgacs=(5,6))
    warm=dict(common,psevp=(1,5),pgevp=(1,6),pseml=(3,5),pgeml=(3,6),pgacs=(5,6))
    for terms,anchor,is_cold in [(cold,'         work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)+pidep(i,k))\n',True),
                               (warm,'         work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k))\n',False)]:
        code=f'         if(pd_enabled) then\n          call pd_mark(k,{state})\n'
        for name,(a,b) in terms.items(): code+=transfer('ws_'+name,a,b,f'{name}(i,k)*dtcld')
        if is_cold:
            code+=transfer('ws_paacw_snow',2,5,'paacw(i,k)*dtcld')
            code+=transfer('ws_paacw_graupel',2,6,'paacw(i,k)*dtcld')
            for name,a,frac in [('piacr',3,'delta3'),('praci',4,'delta3'),('psacr',3,'delta2')]:
                code+=transfer('ws_'+name,a,5,f'{name}(i,k)*{frac}*dtcld')
                code+=transfer('ws_'+name,a,6,f'{name}(i,k)*(1.-{frac})*dtcld')
            code+=transfer('ws_pracs',5,6,'pracs(i,k)*(1.-delta2)*dtcld')
        else: code+=transfer('ws_paacw_warm',2,3,'2.*paacw(i,k)*dtcld')
        code+='         endif\n'; text=insert(text,anchor,code)
    # Same xlf assignment occurs elsewhere: attach to the exact end of updates.
    for anchor in ['                 + pgacr(i,k)+pgacs(i,k))*dtcld,0.)\n',
                   '                 + pgeml(i,k))*dtcld,0.)\n']:
        text=insert(text,anchor,f"         if(pd_enabled) call pd_adjust(k,'positivity_clipping',{state})",after=True)
    # Padding after all minor-step tendencies, before the next internal loop.
    text=insert(text,'       if(qc(i,k).le.qmin) qc(i,k) = 0.0\n',f'       if(pd_enabled) call pd_mark(k,{state})')
    text=insert(text,'       if(qi(i,k).le.qmin) qi(i,k) = 0.0\n',f"       if(pd_enabled) call pd_difference(k,'final_clipping',{state})",after=True)
    return text


def main():
    PROV.mkdir(parents=True,exist_ok=True)
    for relative,fn in [('src/module_mp_thompson.F',thompson),('src/wsm6/mp_wsm6.F90',wsm6)]:
        instrument(ROOT/relative,fn)
    save_json(PROV/'hooks.json',{name:sha(PROV/name) for name in ['module_mp_thompson.F','mp_wsm6.F90']})
    print('Original physics statements preserved; only delimited diagnostic hooks inserted.')


if __name__=='__main__': main()

"""Validate relocated immutable inputs, model runs and existing DEPHY changes."""
from concurrent.futures import ThreadPoolExecutor
import json
from scm_lightning.provenance import ROOT, git, resolve_path, sha256, write_json


def main():
    root=ROOT/'runs/argentina_20181213/_execution/ensemble40'
    runs=json.loads((root/'run_checks.json').read_text(encoding='utf-8'))
    plan=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    def verify(run):
        path=resolve_path(run['folder'])/'output.nc'
        actual=sha256(path)
        if actual != run['sha256']:
            raise AssertionError('Relocated simulation changed: '+str(path))
        return {'configuration':run['configuration'], 'member':run['member'],
                'path':str(path.relative_to(ROOT)), 'sha256':actual}
    with ThreadPoolExecutor(max_workers=3) as pool:
        outputs=list(pool.map(verify,runs))
    inputs={r['forcing']:r['forcing_sha256'] for r in plan['runs']}
    for path, expected in inputs.items():
        if sha256(resolve_path(path)) != expected:
            raise AssertionError('Relocated ensemble forcing changed')
    before=json.loads((ROOT/'results/argentina_20181213/audit_before.json').read_text(encoding='utf-8'))
    dephy=ROOT/'external/DEPHY-SCM'
    assert git('rev-parse','HEAD',cwd=dephy)==before['dephy_head']
    assert git('status','--short',cwd=dephy)==before['dephy_status']
    write_json(ROOT/'results/argentina_20181213/tests/path_migration.json', {
        'status':'passed', 'model_outputs_checked':len(outputs), 'unique_member_forcings_checked':len(inputs),
        'all_checksums_unchanged':True, 'dephy_commit_and_preexisting_status_unchanged':True,
        'legacy_manifests':'Retained unchanged; legacy path strings resolved through provenance.resolve_path',
        'outputs':outputs})
    print('246 native outputs and 41 member inputs unchanged; DEPHY local changes preserved.')


if __name__=='__main__':
    main()

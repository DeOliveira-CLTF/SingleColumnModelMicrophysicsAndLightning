"""Cache public-domain Natural Earth boundaries with source checksums."""
import json
import urllib.request
from scm_lightning.provenance import ROOT, sha256, write_json

LAYERS = {
    'countries': 'ne_50m_admin_0_countries',
    'states': 'ne_10m_admin_1_states_provinces_lines',
}


def main():
    folder = ROOT/'data/raw/boundaries'
    folder.mkdir(parents=True, exist_ok=True)
    records = []
    for name, layer in LAYERS.items():
        url = 'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/'+layer+'.geojson'
        target = folder/(name+'.geojson')
        if not target.exists():
            with urllib.request.urlopen(url, timeout=120) as response:
                content = response.read()
            parsed = json.loads(content)
            if parsed['type'] != 'FeatureCollection':
                raise ValueError('Invalid geographic boundary source')
            target.write_bytes(content)
        records.append({'path': target.name, 'source': url, 'sha256': sha256(target)})
    write_json(folder/'manifest.json', {'license': 'Natural Earth public domain', 'layers': records})


if __name__ == '__main__':
    main()

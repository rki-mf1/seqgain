// Execute the report's own script against a minimal DOM/Plotly test double.
// No browser, network requests, or additional npm dependencies are needed.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const html = fs.readFileSync(process.argv[2], 'utf8');
const payload = html.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const source = html.slice(html.lastIndexOf('<script>') + 8, html.lastIndexOf('</script>'));
const data = JSON.parse(payload), elements = new Map(), plots = new Map();

class Element {
  constructor(tag = 'div') {
    this.tag = tag; this.children = []; this.style = {}; this.textContent = ''; this._value = '';
  }
  set id(id) { this._id = id; elements.set(id, this); }
  get id() { return this._id; }
  set value(value) { this._value = String(value); }
  get value() { return this._value; }
  append(...items) { this.children.push(...items); }
  add(option) { if (!this.children.length) this.value = option.value; this.append(option); }
  replaceChildren(...items) { this.children = items; }
  insertBefore(item, before) { this.children.splice(Math.max(0, this.children.indexOf(before)), 0, item); }
  after() {}
}
const get = id => {
  if (!elements.has(id)) { const e = new Element(); e.id = id; }
  return elements.get(id);
};
get('data').textContent = payload;
get('horizon').value = 2; get('slider').value = 150; get('threshold').value = 10;
const sandbox = {
  document: {getElementById: get, createElement: tag => new Element(tag)},
  Option: function(text, value) { const e = new Element('option'); e.textContent = text; e.value = value; return e; },
  Plotly: {
    newPlot: (id, traces) => plots.set(id, traces),
    react: (id, traces) => plots.set(id, traces),
    relayout: () => {},
  },
};
vm.runInNewContext(source, sandbox, {timeout: 10000});
if (data.bam) {
  const variants = [data.bam, ...(data.bam_filter_variants || [])];
  assert.equal(get('mapq').disabled, variants.length === 1);
  assert.equal(get('duplicates').disabled, variants.length === 1);
  for (const variant of variants) {
    const def = variant.coverage_definition;
    get('mapq').value = def.minimum_mapq;
    get('duplicates').checked = def.include_marked_duplicates;
    get('mapq').onchange();
    const diagnostics = JSON.parse(get('provenance').textContent).reference_coverage;
    assert.deepEqual(diagnostics.definition, def);
    assert.deepEqual(diagnostics.record_counts, variant.record_counts);
    assert.ok(get('warnings').children.some(li => li.textContent.includes(`below MAPQ ${def.minimum_mapq}`)));
    const views = [variant.summary, ...variant.targets];
    const full = t => t.curves.find(r => r.replicate === 0 && r.fraction === 1);
    assert.equal(get('meandepth').textContent, full(views[0]).mean_fragment_depth.toFixed(2) + '×');
    for (let i = 0; i < views.length; i++) {
      get('segment').value = i;
      get('segment').onchange();
      for (const depth of [1, 9, 100]) {
        get('threshold').value = depth; get('threshold').oninput();
        const expected = views[i].curves.filter(r => r.replicate === 0).map(r => r.breadth[depth]);
        assert.deepEqual(Array.from(plots.get('bchart')[0].y), expected);
        const forecast = views[i].model || views[i];
        if (forecast.forecast_breadth[depth]) {
          assert.deepEqual(Array.from(plots.get('bchart').at(-1).y), forecast.forecast_breadth[depth]);
        }
        assert.equal(get('coveragerows').children[i].children[2].textContent,
                     full(views[i]).mean_fragment_depth.toFixed(2)+'×');
        get('slider').value = 100; get('slider').oninput();
        const row = get('rows').children[i];
        assert.equal(row.children[2].textContent, (full(views[i]).breadth[depth]*100).toFixed(1)+'%');
        assert.equal(row.children[2].textContent, row.children[3].textContent);
        assert.ok(get('scenario').textContent.includes(views[i].contig || 'all selected references'));
      }
    }
  }
} else {
  assert.equal(get('bsection').hidden, true);
}
get('horizon').value = 100; get('horizon').oninput();
get('slider').value = 4100; get('slider').oninput();
assert.ok(get('scenario').textContent.includes('41.00×'));
assert.ok(get('scenario').textContent.includes('Long-range scenario'));
if (data.kmers) assert.ok(plots.has('kchart'));
console.log('Report script and filter interactions passed:', process.argv[2]);

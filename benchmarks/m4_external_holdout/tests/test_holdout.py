import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import external_holdout as h


class HoldoutTests(unittest.TestCase):
    def test_exclusion_plan_covers_all_named_main_copies_once(self):
        d = h.Data(h.FROZEN / 'inputs/benchmark')
        plans = h.load_plan('species')
        ids = [x for p in plans for x in p['query_ids']]
        expected = {x for x, r in d.cop.items() if r['culture'] in d.species}
        self.assertEqual(len(plans), 37)
        self.assertEqual(set(ids), expected)
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(sum(p['named_main_queries'] for p in plans), 803)
        for p in plans:
            held = set(p['held_cultures'])
            for c in d.by_species[p['target']]:
                self.assertTrue(d.by_label[d.cul[c]['culture_label']] <= held)

    def test_aliases_are_removed_and_absent_species_is_explicit(self):
        p = {p['target']: p for p in h.load_plan('species')}
        for sp in ('Rhizophagus prolifer', '[Rhizoglomus] vesiculiferum', 'Septoglomus viscosum'):
            self.assertEqual(len(p[sp]['removed_ref_ids']), 2)
        self.assertEqual(p['Ambispora callosa']['removed_ref_ids'], [])
        self.assertEqual(p['Ambispora callosa']['external_status'], 'already_absent_by_label')

    def test_genus_excludes_unnamed_cultures_and_physical_companions(self):
        d = h.Data(h.FROZEN / 'inputs/benchmark')
        plans = h.load_plan('genus')
        for p in plans:
            held = set(p['held_cultures'])
            self.assertTrue({c for c in d.cultures if d.genus[c] == p['target']} <= held)
            self.assertFalse(any(d.genus[c] == p['target'] for c in d.cultures - held))

    def test_pruning_removes_exact_taxa_and_preserves_others(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'tree.newick'
            path.write_text('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);')
            tree = h.pruned_tree(path, {'B', 'D'}, {'A', 'C', 'E'})
            self.assertEqual({n.name for n in tree.get_terminals()}, {'A', 'C', 'E'})
            with self.assertRaises(ValueError):
                h.pruned_tree(path, {'missing'}, {'A', 'B', 'C', 'D', 'E'})

    def test_completed_stage_detects_tampering_and_configuration_changes(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'stage'
            h.stage(path, 'key1', lambda p: (p / 'result.txt').write_text('complete'))
            self.assertTrue(h.valid_stage(path, 'key1'))
            with self.assertRaises(ValueError): h.valid_stage(path, 'key2')
            (path / 'result.txt').write_text('tampered')
            with self.assertRaises(ValueError): h.valid_stage(path, 'key1')

    def test_interrupted_stage_is_restarted_not_accepted(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'stage'; path.mkdir(); (path / 'partial.txt').write_text('incomplete')
            self.assertFalse(h.valid_stage(path, 'key'))
            h.stage(path, 'key', lambda p: (p / 'result.txt').write_text('complete'))
            self.assertFalse((path / 'partial.txt').exists())
            self.assertTrue(h.valid_stage(path, 'key'))


if __name__ == '__main__':
    unittest.main()

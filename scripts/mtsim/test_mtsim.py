"""Unit tests for the mtDNA read simulator's sequence and truth logic."""
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mtsim  # noqa: E402
import score  # noqa: E402

REFERENCE = "ACGTTTTTACGTACGTACNCGTACGTGGGGGACGT"


class EditTests(unittest.TestCase):
    def test_apply_edits_substitution_deletion_insertion(self):
        edited = mtsim.apply_edits("AAACCCGGG", [
            mtsim.Edit(2, 1, "T"), mtsim.Edit(4, 2), mtsim.Edit(8, 0, "AC")])
        self.assertEqual(edited, "ATACGACGG")

    def test_apply_edits_drops_reference_placeholder(self):
        self.assertEqual(mtsim.apply_edits("ACNCT", []), "ACCT")

    def test_overlapping_edits_are_rejected(self):
        with self.assertRaises(ValueError):
            mtsim.apply_edits("AAAAAA", [mtsim.Edit(1, 3), mtsim.Edit(2, 1, "T")])

    def test_parse_cs_long(self):
        edits = mtsim.parse_cs_long("=ACG*tc=AA+gg=T-ca=A", 11)
        self.assertEqual(edits, [mtsim.Edit(14, 1, "C"), mtsim.Edit(17, 0, "GG"),
                                 mtsim.Edit(18, 2)])

    def test_parse_cs_long_skips_placeholder_gap(self):
        self.assertEqual(mtsim.parse_cs_long("=AC-n=GT", 1), [])


class NormalizationTests(unittest.TestCase):
    def test_deletion_moves_left_through_homopolymer(self):
        # Deleting any T of TTTTT gives the same molecule.
        for position in range(4, 9):
            self.assertEqual(mtsim.normalize_indel(REFERENCE, position, 1, ""),
                             ("deletion", 4, "T"))

    def test_insertion_moves_left_through_homopolymer(self):
        for position in range(5, 10):
            self.assertEqual(mtsim.normalize_indel(REFERENCE, position, 0, "T"),
                             ("insertion", 3, "T"))

    def test_repeat_deletion_rotates_allele(self):
        sequence = "GGCATCATCATTT"
        first = mtsim.normalize_indel(sequence, 3, 3, "")
        shifted = mtsim.normalize_indel(sequence, 6, 3, "")
        self.assertEqual(first, shifted)
        self.assertEqual(first, ("deletion", 3, "CAT"))

    def test_normalization_terminates_on_uniform_circle(self):
        self.assertEqual(mtsim.normalize_indel("AAAAAA", 3, 1, "")[0], "deletion")


class TruthTests(unittest.TestCase):
    def test_fraction_sums_over_components_and_skips_numt(self):
        shared = mtsim.Edit(1, 1, "T")
        components = [
            mtsim.Component("a", 0.75, [shared]),
            mtsim.Component("b", 0.25, [shared, mtsim.Edit(3, 1, "A")]),
            mtsim.Component("numt", 0.5, [mtsim.Edit(2, 1, "A")], numt_like=True),
        ]
        records = mtsim.truth_records("ACGTACGT", components)
        by_position = {r["position"]: r for r in records}
        self.assertEqual(by_position[1]["fraction"], 1.0)
        self.assertEqual(by_position[3]["fraction"], 0.25)
        self.assertNotIn(2, by_position)


class ScoringTests(unittest.TestCase):
    def test_a_deletion_with_the_wrong_length_is_not_a_match(self):
        truth = [{"position": 8470, "length": 4977, "fraction": 0.2}]
        exact = [{"left_aligned_start": 8474, "length": 4977, "fraction": 0.2}]
        wrong = [{"left_aligned_start": 8474, "length": 4978, "fraction": 0.2}]
        self.assertEqual(score.score_deletions(exact, truth)["tp"], 1)
        self.assertEqual(score.score_deletions(wrong, truth)["tp"], 0)

    def test_calls_outside_callable_positions_do_not_count_as_false(self):
        truth = {"variants": [{"kind": "snv", "position": 100, "ref": "A", "alt": "G", "fraction": 1.0},
                              {"kind": "snv", "position": 900, "ref": "A", "alt": "G", "fraction": 1.0}],
                 "sample": "t"}
        call = {"status": "PASS", "alt": "G", "fraction": 1.0}
        result = {"precision_calls": {
            "callable_ranges": [[50, 200]],
            "snvs": [{**call, "position": 100}, {**call, "position": 300}],
            "indels": [], "deletions": []}}
        report = score.score(result, truth)["snv"]
        self.assertEqual((report["tp"], report["fp"], report["fn"]), (1, 0, 0))


class ReadSimulationTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(7)
        self.sequence = "".join(rng.choice("ACGT") for _ in range(4000))
        self.cache = {}

    def test_same_seed_gives_same_read(self):
        platform = mtsim.PLATFORMS["ont_r10"]
        template = mtsim.Template(self.sequence, platform, self.cache)
        first = mtsim._simulate_read(template, 10, 800, platform, random.Random(3))
        second = mtsim._simulate_read(template, 10, 800, platform, random.Random(3))
        self.assertEqual(first, second)
        self.assertEqual(len(first[0]), len(first[1]))

    def test_read_can_cross_the_molecule_end(self):
        platform = mtsim.PLATFORMS["hifi"]
        template = mtsim.Template(self.sequence, platform, self.cache)
        read, _ = mtsim._simulate_read(template, 3900, 300, platform, random.Random(1))
        self.assertGreater(len(read), 250)
        self.assertLess(mtsim.hamming_prefix(read, (self.sequence * 2)[3900:]), 40)

    def test_hifi_has_fewer_errors_than_nanopore(self):
        expected = (self.sequence * 2)[0:2000]
        errors = {}
        for name in ("hifi", "ont_r9"):
            platform = mtsim.PLATFORMS[name]
            template = mtsim.Template(self.sequence, platform, {})
            read, _ = mtsim._simulate_read(template, 0, 2000, platform, random.Random(5))
            errors[name] = mtsim.hamming_prefix(read, expected) + abs(len(read) - 2000)
        self.assertLess(errors["hifi"], errors["ont_r9"])

    def test_homopolymer_error_grows_with_length(self):
        platform = mtsim.PLATFORMS["ont_r10"]
        self.assertEqual(mtsim._hp_error(platform, 2), 0.0)
        self.assertLess(mtsim._hp_error(platform, 4), mtsim._hp_error(platform, 8))
        self.assertLessEqual(mtsim._hp_error(platform, 60), platform.hp_cap)


if __name__ == "__main__":
    unittest.main()

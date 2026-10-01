import gzip
import io
import unittest
from fetch_reads import fastq_prefix

READ = b"@read\nGATC\n+\nIIII\n"

class PrefixTests(unittest.TestCase):
    def test_bounded_selection_preserves_exact_bytes_across_chunks(self):
        compressed = gzip.compress(READ * 3)
        output = io.BytesIO()
        result = fastq_prefix((bytes([byte]) for byte in compressed), output, 2, 1024)
        self.assertEqual(result['selected_reads'], 2)
        self.assertEqual(output.getvalue(), READ * 2)

    def test_truncated_record_is_not_a_successful_sample(self):
        with self.assertRaisesRegex(ValueError, 'complete records'):
            fastq_prefix([gzip.compress(READ + b'@partial\n')], io.BytesIO(), 2, 1024)

    def test_malformed_quality_rejected(self):
        with self.assertRaisesRegex(ValueError, 'malformed'):
            fastq_prefix([gzip.compress(b'@read\nGATC\n+\nIII\n')], io.BytesIO(), 1, 1024)

    def test_decompression_budget(self):
        with self.assertRaisesRegex(ValueError, 'byte limit'):
            fastq_prefix([gzip.compress(READ * 100)], io.BytesIO(), 100, 50)

    def test_concatenated_gzip_members(self):
        output = io.BytesIO()
        fastq_prefix([gzip.compress(READ) + gzip.compress(READ)], output, 2, 1024)
        self.assertEqual(output.getvalue(), READ * 2)

if __name__ == '__main__':
    unittest.main()

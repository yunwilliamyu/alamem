use assert_cmd::Command;
use std::fs;
use tempfile::NamedTempFile;

#[test]
fn test_binary_file_output() {
    let temp_file = NamedTempFile::new().expect("Failed to create temp file");
    let output_path = temp_file.path();

    let mut cmd = Command::cargo_bin("alamem")
        .expect("Failed to find binary");

    cmd.arg("-t").arg("1")
       .arg("test_files/NZ_CP013494.1.fna")
       .arg("test_files/hidden_NZ_CP013494.1,2424,15275.fasta")
       .arg(output_path)
       .assert()
       .success();

    let actual_output = fs::read_to_string(output_path)
        .expect("Failed to read the output file");

    let expected_output = fs::read_to_string("test_files/output.txt")
        .expect("Failed to read expected output");

    let actual_normalized = actual_output.replace("\r\n", "\n");
    let expected_normalized = expected_output.replace("\r\n", "\n");

    assert_eq!(actual_normalized, expected_normalized)

}

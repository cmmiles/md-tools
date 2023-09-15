# Builds Rustdocs for deployment to GitHub Pages
cargo doc --no-deps
rm -rf ./docs
echo "<meta http-equiv=\"refresh\" content=\"0; url=md_tools\">" > target/doc/index.html
cp -r target/doc ./docs
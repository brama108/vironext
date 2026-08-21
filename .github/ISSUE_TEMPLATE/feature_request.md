---
name: Feature Request
about: Suggest a new feature or improvement to the ViroNEXT pipeline
title: ''
labels: enhancement
assignees: ''
---

## 🚀 Feature Request — Virus Identification Pipeline

### 📌 Summary
<!-- One–two sentences describing the requested feature or change.
Example: Replace tool X with tool Y for virus genome quality assessment. -->

---

### 🤔 Motivation
- **Problem:** <!-- What issue does the current pipeline have? -->
- **Impact:** <!-- Why is this important (accuracy, runtime, reproducibility, usability)? -->
- **Goal:** <!-- What is the desired outcome of the change? -->

---

### 🛠️ Proposed Solution
- **New Tool/Method:** <!-- Which software/tool/library should be integrated? -->
- **Why This Tool?**
  - Maintained? ✅
  - Widely used in virology? ✅
  - Provides metrics relevant for viral identification? ✅
- **Proposed Workflow Change:**
  - Where in the pipeline would this step occur?
  - Should it *replace* an existing step or be *added in parallel*?

---

### 📂 Implementation Plan
1. <!-- Step 1: e.g., remove tool X from rules/NN_tool.smk -->
2. <!-- Step 2: write new Snakemake rule for tool Y -->
3. <!-- Step 3: update configs (envs/*.yaml, resources, params) -->
4. <!-- Step 4: test on small dataset -->
5. <!-- Step 5: update documentation + example run -->

---

### ✅ Acceptance Criteria
- [ ] Pipeline runs successfully with the new step.
- [ ] Outputs are clearly defined and stored in the results directory.
- [ ] Documentation is updated with usage notes.
- [ ] Benchmarking shows equal or better accuracy/performance than the previous version.

---

### 📚 References
- Tool homepage: <!-- link -->
- Citation: <!-- reference for publication -->
- Related issues/PRs: <!-- cross-links -->

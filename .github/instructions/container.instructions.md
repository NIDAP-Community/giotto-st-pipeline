# File: .github/instructions/container.instructions.md
---
applyTo: "container/**"
---

For container-related changes:
- Keep builds deterministic: pin key versions when feasible; document them.
- Prefer baking R dependencies into the image rather than installing at runtime.
- Ensure `container/run.sh` supports:
  - bind-mounting input/output dirs
  - passing config paths
  - writing logs to the host filesystem
- Any change here must include an updated, copy-pastable build + run example in `container/README.md` and/or `README.md`.

## Distribution (Canonical artifacts)

This repository does NOT store container binaries (Docker images or .sif files) in git.

Canonical distribution uses:
1) A Docker/OCI image published to GitHub Container Registry (GHCR)
2) (Optional) A prebuilt .sif file attached to a GitHub Release for offline / restricted-network HPC use

### 1) Docker image (canonical)
- The canonical runtime artifact is a versioned Docker image in GHCR.
- This is used for local/server Docker execution and as the source for .sif conversion.

### 2) .sif file (optional; release-only)
- For HPC environments where registry pulls are slow or restricted, publish a .sif as a GitHub Release asset.
- Do NOT commit .sif files into the repository.
- Always provide a checksum file (sha256) alongside the .sif.


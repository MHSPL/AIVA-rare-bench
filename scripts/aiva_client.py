"""
Async client for the AIVA variant classification API, with retries,
exponential backoff and rate limit handling.
"""

import asyncio
import sys
from pathlib import Path
from typing import Optional

import httpx

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import (
    AIVA_API_URL, AIVA_API_TOKEN, DEFAULT_ASSEMBLY, CLASSIFICATION_TYPE,
    REQUEST_TIMEOUT, MAX_RETRIES, RETRY_DELAY, RATE_LIMIT_DELAY,
)


class AIVAClient:
    """Async client for AIVA variant classification API."""

    def __init__(
        self,
        base_url: str = AIVA_API_URL,
        token: str = AIVA_API_TOKEN,
        timeout: int = REQUEST_TIMEOUT
    ):
        self.base_url = base_url.rstrip("/")
        self.token = token
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None

    async def __aenter__(self):
        self._client = httpx.AsyncClient(
            base_url=self.base_url,
            headers={
                "Authorization": f"Bearer {self.token}",
                "Content-Type": "application/json"
            },
            timeout=httpx.Timeout(self.timeout)
        )
        return self

    async def __aexit__(self, *args):
        if self._client:
            await self._client.aclose()

    async def classify_variant(
        self,
        variant_data: dict,
        classification_type: str = CLASSIFICATION_TYPE,
        assembly: str = DEFAULT_ASSEMBLY,
        model: Optional[str] = None,
        phenotype_terms: Optional[str] = None,
        description: Optional[str] = None,
        additional_context: Optional[str] = None,
        retries: int = MAX_RETRIES
    ) -> dict:
        """
        Call AIVA API to classify a variant.

        Args:
            variant_data: Dict with variant info (hgvs, rsid, or coordinates)
            classification_type: "acmg" or "amp"
            assembly: "GRCh37" or "GRCh38"
            model: LLM model to use (e.g., "openai/gpt-4o", "google/gemini-2.5-flash")
            phenotype_terms: Disease/phenotype description
            description: Sample description
            additional_context: Additional notes
            retries: Number of retries on failure

        Returns:
            API response dict
        """
        payload = {
            "variant_data": variant_data,
            "classification_type": classification_type,
            "assembly": assembly
        }

        if model:
            payload["model"] = model
        if phenotype_terms:
            payload["phenotype_terms"] = phenotype_terms
        if description:
            payload["description"] = description
        if additional_context:
            payload["additional_context"] = additional_context

        last_error = None
        for attempt in range(retries + 1):
            try:
                response = await self._client.post(
                    "/api/analyze/classify-variant",
                    json=payload
                )
                response.raise_for_status()
                return response.json()

            except httpx.TimeoutException as e:
                last_error = f"Timeout after {self.timeout}s"
                if attempt < retries:
                    delay = RETRY_DELAY * (2 ** attempt)  # Exponential backoff
                    await asyncio.sleep(delay)

            except httpx.HTTPStatusError as e:
                last_error = f"HTTP {e.response.status_code}: {e.response.text[:200]}"

                # Special handling for rate limit (429)
                if e.response.status_code == 429:
                    if attempt < retries:
                        delay = RATE_LIMIT_DELAY * (attempt + 1)  # Increasing delay for rate limits
                        print(f"\n  Rate limited, waiting {delay}s before retry {attempt + 1}/{retries}...")
                        await asyncio.sleep(delay)
                    continue

                # Don't retry on other 4xx errors (client errors)
                if 400 <= e.response.status_code < 500:
                    break

                if attempt < retries:
                    delay = RETRY_DELAY * (2 ** attempt)
                    await asyncio.sleep(delay)

            except Exception as e:
                last_error = str(e)
                if attempt < retries:
                    delay = RETRY_DELAY * (2 ** attempt)
                    await asyncio.sleep(delay)

        return {"error": last_error}

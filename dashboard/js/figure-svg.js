  // ---------- figure export (SVG / TIFF) ----------
  // Each figure card is redrawn as plain SVG primitives from its rendered
  // layout (boxes, circles, text), so the SVG stays editable in vector tools.
  // Exports render on a white background.
  const EXPORT_IGNORE = ".fig-export, .ui-hint, .controls, .seg, select, input, label.field";
  const colorCache = new Map();
  const colorCtx = document.createElement("canvas").getContext("2d", { willReadFrequently: true });

  function rgba(css) {
    if (!css || css === "transparent" || css === "none") return null;
    if (colorCache.has(css)) return colorCache.get(css);
    colorCtx.clearRect(0, 0, 1, 1);
    colorCtx.fillStyle = "#000";
    colorCtx.fillStyle = css;
    colorCtx.fillRect(0, 0, 1, 1);
    const [r, g, b, a] = colorCtx.getImageData(0, 0, 1, 1).data;
    const out = a === 0 ? null : { hex: "#" + [r, g, b].map((v) => v.toString(16).padStart(2, "0")).join(""), a: +(a / 255).toFixed(3) };
    colorCache.set(css, out);
    return out;
  }
  const paint = (attr, c) => (c ? ` ${attr}="${c.hex}"${c.a < 1 ? ` ${attr}-opacity="${c.a}"` : ""}` : ` ${attr}="none"`);
  const xmlEsc = (t) => String(t).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const n2 = (v) => +v.toFixed(2);

  function splitTopLevel(str) {
    const parts = []; let depth = 0, cur = "";
    for (const ch of str) {
      if (ch === "(") depth++;
      if (ch === ")") depth--;
      if (ch === "," && depth === 0) { parts.push(cur.trim()); cur = ""; } else cur += ch;
    }
    if (cur.trim()) parts.push(cur.trim());
    return parts;
  }

  function parseShadows(str) {
    if (!str || str === "none") return [];
    return splitTopLevel(str).map((s) => {
      const color = (s.match(/(rgba?|color|oklab|oklch|lab|lch|hsla?)\([^)]*\)|#[0-9a-f]{3,8}/i) || [])[0];
      const rest = color ? s.replace(color, "") : s;
      const nums = (rest.match(/-?[\d.]+px/g) || []).map(parseFloat);
      return { color: rgba(color), inset: /\binset\b/.test(rest), x: nums[0] || 0, y: nums[1] || 0, blur: nums[2] || 0, spread: nums[3] || 0 };
    });
  }

  function shapeSvg(x, y, w, h, radius, attrs) {
    if (radius === "circle") return `<ellipse cx="${n2(x + w / 2)}" cy="${n2(y + h / 2)}" rx="${n2(w / 2)}" ry="${n2(h / 2)}"${attrs}/>`;
    const r = Math.min(radius, w / 2, h / 2);
    return `<rect x="${n2(x)}" y="${n2(y)}" width="${n2(w)}" height="${n2(h)}"${r > 0 ? ` rx="${n2(r)}"` : ""}${attrs}/>`;
  }

  function boxSvg(cs, r, ox, oy) {
    const out = [];
    const rad = cs.borderTopLeftRadius.trim();
    const radius = rad.endsWith("%") && parseFloat(rad) >= 50 ? "circle" : parseFloat(rad) || 0;
    const x = r.left - ox, y = r.top - oy, w = r.width, h = r.height;
    const shadows = parseShadows(cs.boxShadow).filter((s) => s.blur === 0 && s.spread > 0 && s.color);
    shadows.filter((s) => !s.inset).forEach((s) => {
      const grow = radius === "circle" ? "circle" : radius + s.spread;
      out.push(shapeSvg(x - s.spread + s.x, y - s.spread + s.y, w + 2 * s.spread, h + 2 * s.spread, grow, paint("fill", s.color)));
    });
    const bg = rgba(cs.backgroundColor);
    if (bg) out.push(shapeSvg(x, y, w, h, radius, paint("fill", bg)));
    const grad = cs.backgroundImage && cs.backgroundImage.startsWith("linear-gradient")
      ? (cs.backgroundImage.match(/(rgba?|color|oklab|oklch)\([^)]*\)|#[0-9a-f]{3,8}/gi) || []).map(rgba).filter(Boolean) : [];
    if (grad.length >= 2) {
      const id = "g" + Math.random().toString(36).slice(2, 8);
      out.push(`<defs><linearGradient id="${id}">${grad.map((c, i) => `<stop offset="${n2(i / (grad.length - 1))}" stop-color="${c.hex}"/>`).join("")}</linearGradient></defs>`);
      out.push(shapeSvg(x, y, w, h, radius, ` fill="url(#${id})"`));
    }
    const sides = ["Top", "Right", "Bottom", "Left"].map((sd) => ({
      w: cs["border" + sd + "Style"] === "none" ? 0 : parseFloat(cs["border" + sd + "Width"]) || 0,
      c: rgba(cs["border" + sd + "Color"])
    }));
    const same = sides.every((sd) => sd.w === sides[0].w && sd.c && sides[0].c && sd.c.hex === sides[0].c.hex && sd.c.a === sides[0].c.a);
    if (same && sides[0].w > 0) {
      const bw = sides[0].w;
      out.push(shapeSvg(x + bw / 2, y + bw / 2, w - bw, h - bw, radius === "circle" ? "circle" : Math.max(0, radius - bw / 2), ` fill="none"${paint("stroke", sides[0].c)} stroke-width="${bw}"`));
    } else {
      const [t, rt, b, l] = sides;
      if (t.w && t.c) out.push(shapeSvg(x, y, w, t.w, 0, paint("fill", t.c)));
      if (b.w && b.c) out.push(shapeSvg(x, y + h - b.w, w, b.w, 0, paint("fill", b.c)));
      if (l.w && l.c) out.push(shapeSvg(x, y, l.w, h, 0, paint("fill", l.c)));
      if (rt.w && rt.c) out.push(shapeSvg(x + w - rt.w, y, rt.w, h, 0, paint("fill", rt.c)));
    }
    // Inset shadows sit inside the border, as the browser paints them.
    const [bt, br, bb, bl] = sides.map((sd) => sd.w);
    shadows.filter((s) => s.inset).forEach((s) => {
      const sw = s.spread;
      const ix = x + bl + sw / 2, iy = y + bt + sw / 2, iw = w - bl - br - sw, ih = h - bt - bb - sw;
      if (iw > 0 && ih > 0) out.push(shapeSvg(ix, iy, iw, ih, radius === "circle" ? "circle" : Math.max(0, radius - sw / 2), ` fill="none"${paint("stroke", s.color)} stroke-width="${sw}"`));
    });
    return out.join("");
  }

  const measureCtx = document.createElement("canvas").getContext("2d");
  function fitText(text, font, maxW) {
    measureCtx.font = font;
    if (measureCtx.measureText(text).width <= maxW) return text;
    let lo = 0, hi = text.length;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (measureCtx.measureText(text.slice(0, mid) + "…").width <= maxW) lo = mid; else hi = mid - 1;
    }
    return text.slice(0, lo).trimEnd() + "…";
  }

  function textSvg(node, ox, oy) {
    const raw = node.textContent;
    if (!raw.trim()) return "";
    const el = node.parentElement;
    const cs = getComputedStyle(el);
    const color = rgba(cs.color);
    const fs = parseFloat(cs.fontSize);
    const font = `${cs.fontStyle} ${cs.fontWeight} ${cs.fontSize} ${cs.fontFamily}`;
    // Group words into rendered lines.
    const range = document.createRange();
    const lines = [];
    const addPiece = (text, rect, joined) => {
      const line = lines.find((ln) => Math.abs(ln.top - rect.top) < fs * 0.5);
      if (line) { line.words.push(joined && line.glue ? text : " " + text); line.glue = false; line.right = rect.right; line.bottom = Math.max(line.bottom, rect.bottom); }
      else lines.push({ top: rect.top, bottom: rect.bottom, left: rect.left, right: rect.right, words: [text], glue: false });
    };
    for (const m of raw.matchAll(/\S+/g)) {
      range.setStart(node, m.index);
      range.setEnd(node, m.index + m[0].length);
      const rects = [...range.getClientRects()].filter((q) => q.width > 0);
      if (!rects.length) continue;
      if (rects.every((q) => Math.abs(q.top - rects[0].top) < fs * 0.5)) { addPiece(m[0], rects[0], false); continue; }
      // The browser broke this word across lines (e.g. at a hyphen): split it where it broke.
      let piece = "", pieceRect = null;
      for (let i = 0; i < m[0].length; i++) {
        range.setStart(node, m.index + i);
        range.setEnd(node, m.index + i + 1);
        const r = [...range.getClientRects()].find((q) => q.width > 0);
        if (!r) { piece += m[0][i]; continue; }
        if (pieceRect && Math.abs(r.top - pieceRect.top) >= fs * 0.5) { addPiece(piece, pieceRect, false); piece = ""; pieceRect = null; }
        if (!pieceRect) pieceRect = r;
        piece += m[0][i];
      }
      if (piece && pieceRect) addPiece(piece, pieceRect, false);
    }
    lines.forEach((ln) => { ln.words = [ln.words.join("").trim()]; });
    // Text clipped with an ellipsis on screen is clipped the same way here.
    const clip = el.closest(".lbl, .d");
    const clipRect = clip && getComputedStyle(clip).textOverflow === "ellipsis" && clip.scrollWidth > clip.clientWidth ? clip.getBoundingClientRect() : null;
    return lines.map((ln) => {
      let text = ln.words.join(" ");
      if (cs.textTransform === "uppercase") text = text.toUpperCase();
      else if (cs.textTransform === "lowercase") text = text.toLowerCase();
      if (clipRect) text = fitText(text, font, clipRect.right - ln.left);
      const h = ln.bottom - ln.top;
      const baseline = ln.bottom - (h - fs) / 2 - fs * 0.2;
      return `<text x="${n2(ln.left - ox)}" y="${n2(baseline - oy)}" font-size="${fs}" font-weight="${cs.fontWeight}"${cs.fontStyle !== "normal" ? ` font-style="${cs.fontStyle}"` : ""} font-family="${xmlEsc(cs.fontFamily)}"${cs.letterSpacing !== "normal" && parseFloat(cs.letterSpacing) ? ` letter-spacing="${n2(parseFloat(cs.letterSpacing))}"` : ""}${paint("fill", color)}>${xmlEsc(text)}</text>`;
    }).join("");
  }

  // Inline <svg> charts are copied with their computed paint written as attributes.
  function inlineSvg(svg, ox, oy) {
    const r = svg.getBoundingClientRect();
    const clone = svg.cloneNode(true);
    const src = [...svg.querySelectorAll("*")], dst = [...clone.querySelectorAll("*")];
    src.forEach((el, i) => {
      const cs = getComputedStyle(el), c = dst[i];
      c.removeAttribute("class");
      c.removeAttribute("data-tip");
      // Computed paint is copied below; inline var() styles would not resolve outside the page
      c.removeAttribute("style");
      if (cs.display === "none") { c.remove(); return; }
      const fill = cs.fill.startsWith("url") ? null : rgba(cs.fill);
      c.setAttribute("fill", fill ? fill.hex : "none");
      const fo = parseFloat(cs.fillOpacity) * (fill ? fill.a : 1);
      if (fo < 1) c.setAttribute("fill-opacity", n2(fo));
      const stroke = rgba(cs.stroke);
      if (stroke && parseFloat(cs.strokeWidth) > 0) {
        c.setAttribute("stroke", stroke.hex);
        c.setAttribute("stroke-width", parseFloat(cs.strokeWidth));
        if (stroke.a < 1) c.setAttribute("stroke-opacity", stroke.a);
      }
      if (el.tagName === "text") {
        if (cs.paintOrder && cs.paintOrder !== "normal") c.setAttribute("paint-order", cs.paintOrder);
        if (cs.strokeLinejoin) c.setAttribute("stroke-linejoin", cs.strokeLinejoin);
        c.setAttribute("font-size", parseFloat(cs.fontSize));
        c.setAttribute("font-weight", cs.fontWeight);
        c.setAttribute("font-family", cs.fontFamily);
      }
    });
    return `<g transform="translate(${n2(r.left - ox)},${n2(r.top - oy)})">${clone.innerHTML}</g>`;
  }

  function walkSvg(node, ox, oy, out, root) {
    if (node.nodeType === 3) { out.push(textSvg(node, ox, oy)); return; }
    if (node.nodeType !== 1 || node.matches(EXPORT_IGNORE)) return;
    if (node.tagName.toLowerCase() === "svg") { out.push(inlineSvg(node, ox, oy)); return; }
    const cs = getComputedStyle(node);
    if (cs.display === "none" || cs.visibility === "hidden") return;
    if (node !== root) {
      const r = node.getBoundingClientRect();
      if (r.width > 0 && r.height > 0) out.push(boxSvg(cs, r, ox, oy));
    }
    node.childNodes.forEach((c) => walkSvg(c, ox, oy, out, root));
  }

  function withLightLayout(card, fn) {
    const html = document.documentElement;
    html.classList.add("exporting-root");
    // Let wide tables lay out at full width instead of scrolling.
    const prevWidth = card.style.width;
    const extra = Math.max(0, ...[...card.querySelectorAll(".table-wrap")].map((w) => w.scrollWidth - w.clientWidth));
    if (extra > 0) card.style.width = card.offsetWidth + extra + "px";
    card.classList.add("exporting");
    try { return fn(); } finally {
      card.classList.remove("exporting");
      html.classList.remove("exporting-root");
      card.style.width = prevWidth;
    }
  }

  function figureSvg(card) {
    return withLightLayout(card, () => {
      const pad = 16;
      const r = card.getBoundingClientRect();
      const ox = r.left - pad, oy = r.top - pad;
      const out = [];
      walkSvg(card, ox, oy, out, card);
      const w = Math.ceil(r.width + 2 * pad), h = Math.ceil(r.height + 2 * pad);
      const body = out.join("");
      // Embed the page's fonts (only those the figure uses) so SVG and TIFF keep the same type
      const faces = ($("fontFaces")?.textContent || "").split("\n").filter((r) => body.includes(r.match(/font-family:'([^']+)'/)?.[1] || "\u0000"));
      const defs = faces.length ? `<defs><style>${faces.join("\n")}</style></defs>` : "";
      const svg = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">` +
        `${defs}<rect width="${w}" height="${h}" fill="#ffffff"/>${body}</svg>`;
      return { svg, w, h };
    });
  }


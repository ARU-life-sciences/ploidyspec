// shared definitions for the ploidyspec handbook and its generated species chapter

#let accent = rgb("#2b6cb0")
#let warm = rgb("#c05621")
#let soft = rgb("#ebf4ff")
#let softwarm = rgb("#fffaf0")
#let grey = rgb("#4a5568")
#let green = rgb("#2f855a")

// species cards (filled by handbook_species.typ)
#let confcol(c) = if c == "high" { green } else if c == "medium" { rgb("#b7791f") } else if c == "low" { warm } else { rgb("#a0aec0") }
#let spcard(code, name, known, source, facts, answers, caveats) = block(
  breakable: true, width: 100%, inset: 9pt, radius: 4pt, stroke: 0.5pt + rgb("#cbd5e0"), below: 10pt, [
    #text(weight: "bold", fill: accent, size: 11pt, raw(code)) #h(5pt) #emph(name)
    #if known != "" [ \ #text(size: 8.5pt, fill: grey)[Known: #known#if source != "" [ (#source)]]]
    #v(-0.3em)
    #set text(size: 8.6pt)
    #grid(columns: (3.1cm, 1fr), row-gutter: 4pt, column-gutter: 6pt,
      ..facts.map(f => (text(weight: "bold", f.at(0)), f.at(1))).flatten())
    #v(0.2em)
    #grid(columns: (3.1cm, 1.5cm, 1fr), row-gutter: 5pt, column-gutter: 6pt,
      ..answers.map(a => (
        text(weight: "bold", a.at(0)),
        box(fill: confcol(a.at(1)), inset: (x: 3pt, y: 1.5pt), radius: 2pt, text(size: 7pt, fill: white, weight: "bold", a.at(1))),
        [#a.at(2) #if a.at(3) != "" [\ #text(size: 7.6pt, fill: grey, a.at(3))]],
      )).flatten())
    #if caveats.len() > 0 [
      #v(0.2em)
      #text(size: 8.2pt, fill: warm)[*Caveats:* #caveats.join("; ")]
    ]
  ],
)


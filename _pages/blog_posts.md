---
layout: archive
title: "Blog Posts"
permalink: /blog-posts/
author_profile: true
---

{% include base_path %}
{% for post in site.blog_posts %}
  {% include archive-single.html %}
{% endfor %}